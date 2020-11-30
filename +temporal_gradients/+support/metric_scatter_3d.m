function metric_scatter_3d(coord,z,clim,filename,cmap)
% Wrapper for constructing a 3D scatter plot 
%
%   METRIC_SCATTER_3D(coord,z,clim,filename) constructs
%   a scatter plot N-by-3 matrix coord with intensities z and color limits
%   clim. The plot is stored in file filename.
%
%   METRIC_SCATTER_3D(coord,z,clim,filename,cmap) sets the colormap to
%   N-by-3 matrix cmap. If cmap is not provided then the function defaults
%   to parula.

% Set default colormap.
if ~exist('cmap','var')
    cmap = parula;
end

% Create basic figure.
h.figure = figure('Color','w');
h.axes = axes;
h.sct = scatter3(coord(:,1),coord(:,2),coord(:,3),200,z,'.');
colormap(cmap)
xlabel('Gradient 1')
ylabel('Gradient 2')
zlabel('Gradient 3')

% Set axis properties
set(h.axes                                          , ...
    'PlotBoxAspectRatio'    , [1 1 1]               , ...
    'DataAspectRatio'       , [1 1 1]               , ...
    'XLim'                  , [-4 4]                , ...
    'XTick'                 , [-4 0 4]              , ...
    'YLim'                  , [-4 4]                , ...
    'YTick'                 , [-4 0 4]              , ...
    'ZLim'                  , [-4 4]                , ...
    'ZTick'                 , [-4 0 4]              , ...
    'FontName'              , 'DroidSans'           , ...
    'FontSize'              , 22                    );
caxis(clim)
grid on

% Export figure.
export_fig(char(filename), '-m2', '-png');
close(gcf); 
end