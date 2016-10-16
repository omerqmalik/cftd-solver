function S_fig = plotting_getFigStruct(S_structdata,fig_type,varargin)
    title_str = plotting_getStandardTitle(S_structdata);
    if strcmp(fig_type,'fft_field')
        if length(S_structdata.pump) == 1
            S_fig.type = '2Dfft';
            S_fig.islin = varargin{1};
            S_fig.isdec = varargin{2};

            S_fig.x_label = 'w';
            S_fig.y_label = {['$|FFT[' S_structdata.id '(x_0,t)]|$']};

            title_str = [title_str ', D0=' num2str(S_structdata.pump/S_structdata.th,'%.2f')];
        else
            S_fig.type = '3Dfft';

            S_fig.x_label = 'pump';
            S_fig.y_label = 'w';
            S_fig.z_label = {['$|FFT[' S_structdata.id '(x_0,t)]|$']};
        end
    elseif strcmp(fig_type,'tw')
            S_fig.type = 'tw';
            S_fig.x_label = 't';
            S_fig.y_label = {['$|' S_structdata.id '(x_0,t)|$']};
            if length(S_structdata.pump) == 1            
                title_str = [title_str ', D0=' num2str(S_structdata.pump/S_structdata.th,'%.2f')];
            else
                title_str = [title_str ', D0=' num2str(S_structdata.pump/S_structdata.th,'%.2f') '-' num2str(S_structdata.pump/S_structdata.th,'%.2f')];
            end
    end
    S_fig.title_str = [title_str ' t=' num2str(sum(sum(S_structdata.calc_times))/60/60,'%.2f') 'h' ' (' S_structdata.id S_structdata.type ')'];
    
    S_fig.fname = S_structdata.fname;
end