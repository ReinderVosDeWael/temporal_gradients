function master_figure1()
% Constructs the sub-figures required for figure 1. 
%
%   MASTER_FIGURE1() constructs all the sub-figures required to rebuild
%   figure 1 from Vos de Wael et al., 2020, bioRxiv. Figures are stored
%   inside +temporal_gradients/figures/figure_1/. This function requires
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
    'connectivity_vector_3829', ...
    'surf_lh', ...
    'surf_rh', ...
    'connectivity_distance', ...
    'node_strength', ...
    'c69_20k', ...
    'sc_mask', ...
    'temporalLobe_msk');

% Mask data for visualization
connectivity_vector_3829(~c69_20k.mask) = -inf;

% Basic connectivity figures
figure_dir = package_dir + fs + 'figures' + fs + 'figure_1' + fs;
if ~exist(figure_dir, 'dir')
    mkdir(figure_dir)
end

% Builds figures - consult nested figure functions for details. 
build_basic_connectivity_surface(connectivity_vector_3829,surf_lh,surf_rh, figure_dir);
build_affinity_matrix(sc_mask,temporalLobe_msk,c69_20k.mask,figure_dir);
build_scree_plot(gm_hcp_discovery.lambda{1}, figure_dir + 'left_scree.png');
build_gradient_surfaces([gm_hcp_discovery.aligned{1};gm_hcp_discovery.aligned{2}],surf_lh,surf_rh,temporalLobe_msk, figure_dir);
build_graph_scatters(gm_hcp_discovery,{connectivity_distance.hcp_discovery,node_strength.hcp_discovery}, ...
    ["Connectivity Distance","Node Strength"], figure_dir);
build_graph_scatters_3d(gm_hcp_discovery,{connectivity_distance.hcp_discovery,node_strength.hcp_discovery}, ...
    [10,20; 0.5e7,2e7],{parula,parula},["Connectivity Distance","Node Strength"], figure_dir);
end

%% Figure builders
function build_basic_connectivity_surface(vec, surf_lh, surf_rh, figure_dir)
% Builds surfaces of the connectivity of the left temporal pole.
%
%   BUILD_BASIC_CONNECTIVITY_SURFACE(vec, surf_lh, surf_rh, figure_dir)
%   plots data vector vec (here, left temporal pole connectivity) onto the
%   left and right hemispheric surfaces surf_lh and surf_rh. The figure is
%   save as figure_dir/leftPoleConnectivity.png

% Use BrainSpace to construct a surface plot. 
obj = plot_hemispheres(vec,{surf_lh,surf_rh});

% Chang the color limits and map.
set(obj.handles.axes,'CLim',[-3.5 4.5]);
obj.handles.cb.Ticks = [-3.5 4.5];
obj.handles.cb.Limits = [-3.5 4.5];
cmap = hot(256);
cmap = cmap(1:220,:);
colormap([cmap;.7 .7 .7]);

% Add a sphere in the seed. 
coord = surf_lh.coord(:,3829);
for ii = 1:2
    axes(obj.handles.axes(ii)); hold on 
    scatter3(obj.handles.axes(ii),coord(1),coord(2),coord(3),2000,[0 0 0],'.')
end

% Export figure. 
export_fig(char(figure_dir + 'leftPoleConnectivity.png'),'-m2','-png');
close(gcf);
end


function build_affinity_matrix(sc,mask_temporal,mask_midline,figure_dir)
% Display the affinity matrix.
%
%   BUILD_AFFINITY_MATRIX(sc, mask_temporal, mask_midline, figure_dir)
%   plots the affinity matrix of the temporal lobe connectivity (sc) to the
%   whole brain. mask_temporal contains trues for the temporal lobe and
%   mask_midline contains trues for everything but the midline. The output
%   figure is stored as figure_dir/affinity.png. 
%
%   As BrainSpace does not provide an option to return the computed
%   affinity matrix, the affinity matrix (sparsity=75, kernel='cs') is
%   built from scratch here. 

% Left hemispheric connectivity
sc_mask = sc(mask_midline,mask_temporal);
sc_left = sc_mask(1:end/2,1:end/2);

% Sparse connectivity
sc_left(sc_left < prctile(sc_left,75)) = 0;

% Affinity matrix
cosine_similarity = 1-squareform(pdist(sc_left','cosine'));

% Build the figure
h.fig = figure('color','w');
h.ax = axes();
h.img = imagesc(cosine_similarity);
h.cb = colorbar;
set(h.ax                                , ...
    'Box'               , 'off'         , ...
    'DataAspectRatio'   ,  [1 1 1]      , ...
    'PlotBoxAspectRatio', [1 1 1]       , ...
    'Xtick'             , []            , ...
    'YTick'             , []            , ...
    'Visible'           , 'off'         ); 
set(h.cb                                , ...
    'Position'          , [ .85 .3 .03, .4], ...
    'Ticks'             , [0 1]         , ...
    'FontName'          , 'DroidSans'   , ...
    'FontSize'          , 36            );
export_fig(char(figure_dir + 'affinity.png'),'-m2','-png');

end


function build_scree_plot(lambda,name)
% Builds a scree plot
%
%   BUILD_SCREE_PLOT(lambda,name) builds a scree plot from eigenvalue
%   vector lambda. The figure is stored with filename name. 

% Build scree plot using BrainSpace.
h = scree_plot(lambda);

% Adjust axes ticks, limits, and font. 
set(h.axes,'XTick',[0 40],'YTick',[0 .35], 'FontName', 'DroidSans',  ...
    'XLim', [0, 40], 'YLim', [0, .35], 'FontSize', 38)

% Increase marker size.
set(h.plot,'Marker','.','MarkerSize',25);

% Export figure. 
export_fig(char(name),'-m2','-png');
close(gcf);

end


function build_gradient_surfaces(gradients,surf_lh,surf_rh,temporalLobe_msk,figure_dir)
% Builds the surfaces displays of the first three gradients.
%
%   BUILD_GRADIENT_SURFACES(gradients,surf_lh,surf_rh,temporalLobe_msk,
%   figure_dir) displays the first three column vectors of gradients on the
%   surfaces surf_lh and surf_rh. A temporal lobe mask (temporalLobe_msk)
%   is required to plot the data at the correct location. The figure is
%   stored in figure_dir/gradients.png

% Create a temporal lobe 'parcellation'.
fake_parcellation = zeros(20000,1);
fake_parcellation(temporalLobe_msk) = 1:3428;

% Use BrainSpace to plot the hemispheres. 
obj = plot_hemispheres(gradients(:,1:3),{surf_lh,surf_rh},'LabelText',"Gradient " + (1:3), ...
    'Parcellation',fake_parcellation);

% Adjust colormap and color limits. 
colormap([.7 .7 .7; parula])
lim = [-3.6 3.6; -2.4 2.1; -.7 1.4];
obj.colorlimits(lim);

% Export figure. 
export_fig(char(figure_dir + 'gradients.png'),'-m2','-png');
close(gcf);
end


function build_graph_scatters(GM,data,modality_name, figure_dir)
% Build connectivity distance and degree centrality 2D scatter plots.
%
%   BUILD_GRAPH_SCATTERS(GM,data,modality_name,figure_dir) builds scatter
%   plots of eccentricity versus vectors in data. Scatter plots are stored
%   in figure_dir/{side}_{modality}.png where side is left/right and
%   modality is connectivity distance or node strength. 

hemi_name = ["left", "right"];
import temporal_gradients.support.eccentricity
import temporal_gradients.support.metric_scatter

for hemi = 1:2
    % Get indices for the left/right hemispheric data.
    if hemi == 1
        idx = 1:1714;
    else
        idx = 1715:3428;
    end
        for modality = 1:2
            % Have the y label and limits depend on the modality.
            if modality == 1
                ylim = [10 20];
                ylab = 'Connectivity Distance';
            else
                ylim = [0 3.5]*10e6;
                ylab = 'Node Strength';
            end
            % Get the file name.
            name = figure_dir + hemi_name(hemi) + "_" +  ...
                modality_name(modality) + ".png";
            % Plot data and export. 
            metric_scatter(eccentricity(GM.aligned{hemi}(:,1:3)), ...
                           data{modality}(idx), ...
                           [0,5],ylim,'Eccentricity',ylab,true,name);
        end
end
end


function build_graph_scatters_3d(GM,data,clim,cmap,modality_name, figure_dir)
% Build connectivity distance and degree centrality 3D scatter plots.
%
%   BUILD_GRAPH_SCATTERS_3D(GM,data,clim,cmap,modality_name,figure_dir)
%   builds scatter plots of vectors in data in the manifold space described
%   by the first three gradients in the GradientMaps object GM. Variables
%   clim and cmap provide the color limits and colormap, respectively.
%   Scatter plots are stored in figure_dir/{side}_{modality}.png where side
%   is left/right and modality is connectivity distance or node strength.

hemi_name = ["lh", "rh"];
import temporal_gradients.support.eccentricity
import temporal_gradients.support.metric_scatter_3d

for hemi = 1:2
    % Get indicies for left/right hemisphere. 
    if hemi == 1
        idx = 1:1714;
    else
        idx = 1715:3428;
    end
    for modality = 1:numel(data)
        % Get filename.
        name = figure_dir + "scatter3d_" + ...
            hemi_name(hemi) + "_" + modality_name(modality) +".png";
        % Plot and export figure.
        metric_scatter_3d(GM.aligned{hemi}, data{modality}(idx), ...
            clim(modality,:),replace(name,'T1w/T2w.png','T1wT2w.png'), cmap{modality});
    end
end
end