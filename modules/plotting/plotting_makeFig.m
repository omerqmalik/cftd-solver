function fig = plotting_makeFig
    fig = figure;
    pos = fig.Position;
    pos = [pos(1) pos(2) 2*pos(3) 2*pos(4)];
    fig.Position = pos;
end