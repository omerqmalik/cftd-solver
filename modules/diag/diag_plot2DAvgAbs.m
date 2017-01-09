function f = diag_plot2DAvgAbs(S_structdata,S_fig)
    if strcmp(S_fig.type,'DmnAVGABS')
        f = plotting_plotMatrix(S_structdata.Y,S_fig);
    end
end