function S_fig = plotting_getFigStruct(S_structdata,fig_type,opt_islin,opt_isdec)  %This function will always return a 1-element struct S_fig
    title_str  = plotting_getStandardTitle(S_structdata);
    S_fig.type = fig_type;
    
    if strcmp(fig_type,'3DfieldFFT')
        S_fig.x_label = 'pump';
        S_fig.y_label = 'w';
        S_fig.z_label = {['$|FFT[' S_structdata.id '(x_0,t)]|$']};
        
        title_str = [title_str ', x0=' num2str(S_structdata.x0)];
    elseif strcmp(fig_type,'3DcoeffsFFT')
        S_fig.x_label = 'pump';
        S_fig.y_label = 'w';
        S_fig.z_label = {['$|FFT[\sum ' S_structdata.id '(t)]|$']};
    elseif strcmp(fig_type,'2DfieldFFT')
        S_fig.islin = opt_islin;
        S_fig.isdec = opt_isdec;

        S_fig.x_label = 'w';
        S_fig.y_label = {['$|FFT[' S_structdata.id '(x_0,t)]|$']};

        title_str = [title_str ', D0=' num2str(S_structdata.pump/S_structdata.th,'%.2f') '(pstep ' num2str(S_structdata.psteps) '), x0=' num2str(S_structdata.x0)];
    elseif strcmp(fig_type,'2DcoeffsFFT')
        S_fig.islin = opt_islin;
        S_fig.isdec = opt_isdec;
        S_fig.x_label = 'w';
        
        if S_fig.islin == 1
            S_fig.colormap = 'lines';
        elseif S_fig.islin == 0
            S_fig.colormap = 'autumn';
        end

        lgnd = cell(1,S_structdata.nCF);
        for i = 1:S_structdata.nCF
            lgnd{i} = ['$|FFT[' S_structdata.id '_{' num2str(i) '}(t)]|$'];
        end
        S_fig.lgnd = lgnd;
        
        title_str = [title_str ', D0=' num2str(S_structdata.pump/S_structdata.th,'%.2f')];
    elseif strcmp(fig_type,'DcoeffsFFT')
        S_fig.islin = opt_islin;
        S_fig.isdec = opt_isdec;
        S_fig.x_label = 'w';
        
        if S_fig.islin == 1
            S_fig.colormap = 'lines';
        elseif S_fig.islin == 0
            S_fig.colormap = 'autumn';
        end

        G = S_structdata.Dsave - 2*S_structdata.nCF;
        lgnd = cell(1,length(G));
        for i = 1:length(G)
            [a,b] = ind2sub([S_structdata.nCF S_structdata.nCF],G(i));
            lgnd{i} = ['$' S_structdata.id '_{' num2str(b) ',' num2str(a) '}(t)]|$'];   %NOTE: b,a are row/column flipped because the D matrix is naturally transposed
        end
        S_fig.lgnd = lgnd;
        
        title_str = [title_str ', D0=' num2str(S_structdata.pump/S_structdata.th,'%.2f')];
    elseif strcmp(fig_type,'2DfieldTW')
        S_fig.islin = 1;
        S_fig.x_label = 't';
        S_fig.y_label = {['$|' S_structdata.id '(x_0,t)|$']};
        if length(S_structdata.pump) == 1            
            title_str = [title_str ', D0=' num2str(S_structdata.pump/S_structdata.th,'%.2f')];
        else
            title_str = [title_str ', D0=' num2str(S_structdata.pump/S_structdata.th,'%.2f') '-' num2str(S_structdata.pump/S_structdata.th,'%.2f') ', x0=' num2str(S_structdata.x0)];
        end
    elseif strcmp(fig_type,'DmnAVGABS')
        S_fig.x_label = 'n';
        S_fig.y_label = 'm';
        S_fig.z_label = {'$\bar{|D_mn|}$'};
        
        title_str = [title_str ', D0=' num2str(S_structdata.pump/S_structdata.th,'%.2f')];
    elseif strcmp(fig_type,'2DcoeffsSPCFLD')
        S_fig.islin = 1;
        S_fig.isdec = 0;
        
        S_fig.x_label = '$x$';
        S_fig.y_label = '$E(x,t_0)$';
        
        S_fig.colormap = 'jet';
        
        if length(S_structdata.t_ind) == 1
            title_str = [title_str ', D0=' num2str(S_structdata.pump/S_structdata.th,'%.2f') ', tstep=' num2str(S_structdata.t(S_structdata.t_ind)) ', x0=' num2str(S_structdata.x0)];
        else
            title_str = [title_str ', D0=' num2str(S_structdata.pump/S_structdata.th,'%.2f') ', tstep=' num2str(S_structdata.t(S_structdata.t_ind(1))) '-' num2str(S_structdata.t(S_structdata.t_ind(end))) ', x0=' num2str(S_structdata.x0)];
        end
    end
    
    if size(S_structdata.calc_times,2) > 1      %multiple pump steps
        S_fig.title_str = [title_str ' t=' num2str(sum(sum(S_structdata.calc_times))/60/60,'%.2f') 'h' ' (' S_structdata.basis_type ',' S_structdata.id S_structdata.type ')'];
        fname = [fig_type '_' S_structdata.id S_structdata.type ... 
            '_p' num2str(S_structdata.pump(1)/S_structdata.th) '-' num2str(S_structdata.pump(end)/S_structdata.th) ...
            '_ps' num2str(S_structdata.psteps(1)) '-' num2str(S_structdata.psteps(end)) ...
            '_t' num2str(S_structdata.t(1)) '-' num2str(S_structdata.t(end))];
    else    %single pump step
        S_fig.title_str = [title_str ' t=' num2str(sum(S_structdata.calc_times)/60,'%.2f') 'm' ' (' S_structdata.basis_type ',' S_structdata.id S_structdata.type ')'];
        fname = [fig_type '_' S_structdata.id S_structdata.type ...
            '_ps' num2str(S_structdata.psteps) ...
            '_p' num2str(S_structdata.pump(1)/S_structdata.th) ...
            '_t' num2str(S_structdata.t(1)) '-' num2str(S_structdata.t(end))];
    end
    
    append_str = '';
    if nargin > 2
        if ~opt_islin
            append_str = [append_str '_log'];
        end
        if opt_isdec
            append_str = [append_str '_d'];
        end
    end
    S_fig.fname = [fname append_str];
end