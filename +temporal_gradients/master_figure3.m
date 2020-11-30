function master_figure3()
% Constructs the sub-figures required for figure 3. 
%
%   MASTER_FIGURE3() constructs all the sub-figures required to rebuild
%   figure 3 from Vos de Wael et al., 2020, bioRxiv. Figures are stored
%   inside +temporal_gradients/figures/figure_3/. This function requires
%   that the data was already downloaded (see
%   temporal_gradients.download_data). 
%
%   Consult the documentation of the nested functions for details on each
%   sub-figure.
%
%   Written by Reinder Vos de Wael, MICA lab, Nov 2020
%   For further details see our <a
%   href="https://github.com/MICA-MNI/micaopen/tree/master/temporal_gradients">Github</a> page.  

% Load data
fs = string(filesep);
package_dir = regexp(mfilename('fullpath'),'.*\+temporal_gradients','match','once');
load(package_dir + fs + "data" + fs + "figure_data.mat", ...
    'gm_hcp_discovery', ...
    'surf_lh', ...
    'surf_rh', ...
    'temporalLobe_msk', ...
    'yeo', ...
    'yeo_tl_predicted', ...
    'kfold_fc_r', ...
    'mics_fc_r', ...
    'repl_fc_r', ...
    'include', ...
    'r', ...
    'r_ho');

% Build figures
figure_dir = char(package_dir + fs + 'figures' + fs + 'figure_3' + fs);
if ~exist(figure_dir, 'dir')
    mkdir(figure_dir)
end
yeo_figures(gm_hcp_discovery, ...
    temporalLobe_msk(1:end/2), ...
    include.l, ...
    double(yeo(1:end/2)), ...
    yeo_tl_predicted{1}, ...
    surf_lh, ...
    figure_dir);
decision_tree_figures({kfold_fc_r,repl_fc_r,mics_fc_r}, ...
    cat(3,r,r_ho), ...
    {surf_lh,surf_rh}, ...
    temporalLobe_msk, ...
    figure_dir);
end

function yeo_figures(gm, mask, include, canonical, predicted, surface, figure_dir)
% Plots all figures associated with the yeo sub-figure.
%
%   yeo_figures(gm, mask, include, canonical, predicted, surface,
%   figure_dir) plots canonical and predicted Yeo networks on the surface,
%   as well as a 3D scatter plot of the canonical Yeo network in 3D
%   manifold space. GM is the GradientMaps object, mask is a temporal_lobe
%   mask, include is a vector containing the overlap between our temporal
%   lobe mask and the vertices included in the Yeo networks (there's some
%   mismatch at the midline), canonical are the canonical Yeo networks,
%   predicted are the predicted Yeo networks, surface is a cortical
%   surface, and figure_dir is the directory where the figures will be stored.

% Define the colormap.
yeo_cmap = [178.5 178.5 178.5; 
                120 18 134;
                70 130 180;
                0 118 14; 
                196 58 250;
                220 248 164; 
                230 148 34;
                205 62 78]/255;

% Create a temopral lobe 'parcellation'.
fake_parcellation = zeros(10000,1);
idx = find(mask);
fake_parcellation(idx(include)) = 1:sum(include);

%% Plot the canonical networks.
obj = plot_hemispheres(canonical,surface, ...
    'labeltext','Yeo Networks');
obj.colormaps(yeo_cmap);
obj.colorlimits([0,7])
export_fig([figure_dir, 'yeo_canonical.png'], '-png', '-m2');
close(gcf);

%% Plot the predicted networks.
obj = plot_hemispheres(predicted,surface,'parcellation',fake_parcellation, ...
    'labeltext',{{'Predicted', 'Networks'}});
obj.colormaps(yeo_cmap);
obj.colorlimits([0,7])
export_fig([figure_dir, 'yeo_predicted.png'], '-png', '-m2');
close(gcf);

%% Plot the 3D scatter plot. 
import temporal_gradients.support.metric_scatter_3d
metric_scatter_3d(gm.aligned{1}(include,1:3),canonical(idx(include)), [0,7], ...
    [figure_dir, '/yeo_scatter3.png'], yeo_cmap);
end


function decision_tree_figures(subject_level,vertex_level,surfaces, temporalLobe_msk, figure_dir)
% Plots all figures associated with the decision tree sub-figure.
%
%   decision_tree_figures(subject_level,vertex_level,surfaces,
%   temporalLobe_msk, figure_dir) plots histograms of the subject_level
%   correlations of all three datasets across both hemispheres. Next it
%   plots the vertex-level correlations on the cortical surfaces. A
%   temporal lobe mask must be provided. Resulting figures are stored in
%   figure_dir.

%% Subject-level correlations.
h.fig = figure('Color','w','Units','Normalized','Position',[0 0 .7 .7]);

% Plot a histogram for each hemisphere/site.
for ii = 1:3
    for hemi = 1:2        
        idx = (hemi-1)*3+ii;
        h.ax(ii,hemi) = subplot(2,3,idx);
        h.hist(ii,hemi) = histogram(subject_level{ii}(:,hemi),10);
        xlabel('Correlation (r)');
    end
end

% Set figure properties.
set(h.ax                                    , ...
    'Box'                   , 'off'         , ...
    'PlotBoxAspectRatio'    , [1,1,1]       , ...
    'FontSize'              , 16            , ...
    'FontName'              , 'DroidSans'   , ...
    'XTick'                 , [0.25, 0.65]   , ...
    'XLim'                  , [0.25, 0.65]   , ...
    'YTick'                 , [0, 20]       , ...
    'YLim'                  , [0, 20]       );
set(h.hist                                  , ...
    'FaceColor'             , [.7 .7 .7]    );
for ii = 1:3
    h.ax(ii,2).Position = h.ax(ii,1).Position + [0 -.5 0 0];
end

% Export
export_fig([figure_dir, 'subjectwise_correlations.png'],'-png','-m2')
close(gcf);

%% Vertexwise correlations

% Plot vertex level correlations onto the surface.
fake_parcellation = zeros(20000,1);
fake_parcellation(temporalLobe_msk) = 1:3428;
plot = reshape(vertex_level,3428,3);
plot(plot<0) = 0.0001; % A bit above 0 otherwise gray will be assigned to these vertices too. 
obj = plot_hemispheres(plot, surfaces, 'parcellation',fake_parcellation, ...
    'labeltext', {{'HCP', 'Discovery'},{'HCP','Replication'},'MICS'});

% Set colormap and color limits. 
obj.colormaps([.7 .7 .7; parula(10000)]);
obj.colorlimits([0 .9])

% Export
export_fig([figure_dir, 'vertexwise_correlations.png'], '-png', '-m2');
close(gcf);
end
