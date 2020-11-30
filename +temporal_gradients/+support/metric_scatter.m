function metric_scatter(x,y,xlim,ylim,xlab,ylab,add_corr,filename)
% Wrapper for constructing a 2D scatter plot 
%
%   METRIC_SCATTER(x,y,xlim,ylim,xlab,ylab,add_cor,filename) constructs a
%   2D scatter plot of x and y with limits xlim and ylim, and labels xlab
%   and ylab. If add_corr is true, then a Pearson correlation value is
%   added to the plot. The plot is stored in file filename.

% Construct the basic figure.
h.figure = figure('Color','w');
h.axes = axes;
h.sct = scatter(x,y,'.');

% Add a correlation value.
if add_corr
    r = corr(x(:),y(:),'rows','pairwise'); 
    h.txt = text(0.01,0.9,sprintf('r = %0.2f',r), ...
        'Units','normalized', 'FontName', 'DroidSans', ...
        'FontSize', 14);
else
    h.txt = gobjects();
end

% Set the axis properties.
set(h.axes                          , ...
    'PlotBoxAspectRatio', [2 1 1]   , ...
    'XLim'              , xlim       , ...
    'YLim'              , ylim      , ...
    'XTick'             , xlim      , ...
    'YTick'             , ylim      , ...
    'FontName'          , 'DroidSans', ...
    'FontSize'          , 14         );

% Set the marker properties.
set(h.sct, 'SizeData', 20, 'MarkerEdgeColor', [.1 .1 .1]);

% Add x and y labels.
h.txt(end+1) = text(0.5, -0.1, xlab, 'HorizontalAlignment','center', 'FontName', 'DroidSans', ...
    'FontSize', 14, 'Units', 'Normalized');
h.txt(end+1) = text(-0.08, 0.5, ylab, 'HorizontalAlignment','center', 'FontName', 'DroidSans', ...
    'FontSize', 14, 'Units', 'Normalized', 'Rotation', 90);

% Export the figure. 
export_fig(char(filename), '-m2', '-png');
close(gcf); 
end