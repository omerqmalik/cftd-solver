function [d0,mmef,fnef,fnfval] = salt_calcNLSolution(S_fine,S_coredata,S_saltparameters)
    
    d0min = S_fine.th(1)*1.001;                %minimum test/seed d0, x1.001 usually, x1.03 for ECHO_EBONY
    d0max = S_coredata.prel_f*d0min;           %maximum test/seed d0
    dd0 = (d0max-d0min)/S_coredata.pump_sz;    %step size in d0
    d0 = d0min:dd0:d0max;                      %vector of all test/seed d0 values
    
    mmef = zeros(1,length(d0));
    fnef = zeros(length(S_fine.k),length(d0));
    fnfval = zeros(length(S_fine.k),length(d0));

    i = 1;
    pump_is_on = 1;
    no_converge = 0;
    while pump_is_on        %loop through all seed d0 values
        
        %Select appropriate pump step
        if no_converge == 1
            d0(i) = [];
            if i > length(d0)
                break;
            end
        end
        if i == 1
            d0_this = d0(1);
            on_modes = 1;
        else
            is_mode_on = (S_fine.th < d0(i));
            if (sum(is_mode_on) - on_modes >= 2)
                d0 = [d0(1:i-1) (S_fine.th(on_modes+1)+S_fine.th(on_modes+2))/2 d0(i:end)];
                d0_this = d0(i);
                on_modes = sum(S_fine.th < d0(i));
            elseif (sum(is_mode_on) - on_modes < 0)
                no_converge = 1;
                continue;
            else
                d0_this = d0(i);
                on_modes = sum(is_mode_on);
            end
            
            if i == length(d0)
                pump_is_on = 0;
            end
        end

        % Determine how many modes we need for this pump step
        S_crit.modes      = (S_fine.th < d0_this);
        S_crit.k          = S_fine.k(S_crit.modes);
        S_crit.maxpolepos = S_fine.maxpolepos(S_crit.modes);
        S_crit.ord        = S_fine.ord(S_crit.modes);
        S_crit.a          = S_fine.tves(:,S_crit.modes)*S_saltparameters.dtves;     %"direction" vectors for modes < critical value 
        
        %Make sure a-vectors are in the right column
        if i~=1
            newest_mode = 1:length(S_crit.k);
            for aitr = 1:size(S_nlsol.a,2)   %For each mode that has turned ON
                [~,modePos] = min(abs(S_crit.k-S_nlsol.k(aitr)));  %Of all modes turned ON, select one closest to S_nlsol.k
                S_crit.a(:,modePos) = S_nlsol.a(:,aitr);   %The closest mode to S_nlsol.k is best defined by the new direction vector
                newest_mode(newest_mode == modePos) = [];
                if max(S_crit.a(:,modePos))<1e-3
                    S_crit.a(:,modePos)= S_crit.a(:,modePos)/max(abs(S_crit.a(:,modePos)))*0.005;       %Normalization?
                end
            end
        end
        clear('atmp','aSize','aitr','modePos','aLen')

        g = sprintf('%f ', d0);
        fprintf('Current pump vector: %s\n',g);
        fprintf('Current pump step: (%d, %f)\n',i,d0_this);
        
        no_converge = 0;
        while 1
            fprintf('Running nonlinear solver\n');
            g = sprintf('%f ', S_crit.k);
            fprintf('Current ON frequencies: (%d, %s)\n',length(S_crit.k),g);
            
            [S_nlsol,holeburning,mmef(i)] = salt_callNLFsolve(d0_this,on_modes,S_crit,S_coredata,S_saltparameters);       %Calculate the nonlinear term, send only info of modes below critical value (those have turned on)
            
            if mmef(i) < 1
                if mmef(i) == -3       %Trust region radius too small
                    fprintf('Did not converge because of trust issues. Skip this pump step.\n');
                    no_converge = 1;    %this is a special case, flag it down
                    break;
                else
                    fprintf('Did not converge because too many modes turned on. Turn newest mode off and try again.\n');
                    S_crit.a(:,newest_mode) = [];
                    S_crit.k(newest_mode)   = [];
                    S_crit.maxpolepos(newest_mode) = [];
                    S_crit.ord(newest_mode) = [];
                    on_modes = on_modes - 1;
                end
                fprintf('Run nonlinear solver again \n');
            else
                fprintf('Nonlinear solver successful!\n\n\n');
                break;
            end
        end
        
        if no_converge == 1         %Did we have trust region radius issues?
            continue;               %Then skip to the next d0
        end;
            
        fn = [S_coredata.data_dir '/nlsolver_' num2str(i) '.mat'];
        S_nlsol.d0_this = d0_this;
        save(fn,'S_nlsol');

        %----------------------- TH2 finder -------------------------------
        fprintf('Run threshold solver.\n\n');
        S_fine = salt_refineThresholds(S_fine,S_coredata,S_saltparameters,holeburning);
                
        fn = [S_coredata.data_dir '/thsolver_' num2str(i) '.mat'];
        save(fn,'S_fine');
        
        clear('holeburning');
        i = i + 1;
    end
    save([S_coredata.data_dir '/pump.mat'],'d0');
end