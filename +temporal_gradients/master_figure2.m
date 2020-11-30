function master_figure2()
% Constructs the sub-figures required for figure 2. 
%
%   MASTER_FIGURE2() constructs all the sub-figures required to rebuild
%   figure 2 from Vos de Wael et al., 2020, bioRxiv. Figures are stored
%   inside +temporal_gradients/figures/figure_2/. This function requires
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
import temporal_gradients.support.eccentricity
fs = string(filesep);
package_dir = regexp(mfilename('fullpath'),'.*\+temporal_gradients','match','once');
load(package_dir + fs + "data" + fs + "figure_data.mat", ...
    'gm_hcp_discovery', ...
    'surf_lh', ...
    'surf_rh', ...
    'temporalLobe_msk', ...
    'microstructural_features', ...
    'bigbrain', ...
    'sjh');

% Do some data formatting. 
microstructure.lh = table2array(microstructural_features.hcp_discovery(:,1:3));
microstructure.rh = table2array(microstructural_features.hcp_discovery(:,4:6));
hemi = ["lh","rh"];
for ii = 1:2
    h = hemi{ii};
    micro_all.(h) = [microstructure.(h), bigbrain(:,ii)];
end

% Build figures
figure_dir = char(package_dir + fs + 'figures' + fs + 'figure_2' + fs);
if ~exist(figure_dir, 'dir')
    mkdir(figure_dir)
end

% Preset color limits and maps. 
clims = [-0.15, 0.15; 2, 3.5; 1.2, 2];
cmaps = {[.7 .7 .7; parula(256)],[.7 .7 .7; parula(256)],[.7 .7 .7; parula(256)]};

% Build figures. 
build_graph_scatters(gm_hcp_discovery, ...
    microstructure, ...
    {'curvature','thickness','t1wt2w'}, ...
    figure_dir);
build_graph_scatters_3d(gm_hcp_discovery, ...
    micro_all, ...
    [clims; -0.15,0.15], ...
    [cmaps,[.7 .7 .7; greensequential(256)]], ...
    {'Curvature','Thickness','T1w/T2w','MPC'}, ...
    figure_dir);
build_graph_bigbrain(gm_hcp_discovery, ...
    bigbrain, ...
    sjh, ...
    figure_dir);
build_surfaces_tl([[microstructure.lh; microstructure.rh], bigbrain(:)], ...
    {surf_lh,surf_rh}, ...
    [clims;-0.15,0.15], ...
    [cmaps,{[.7 .7 .7; greensequential(256)]}], ...
    {'Curvature','Thickness','T1w/T2w','MPC'}, ...
    temporalLobe_msk, ...
    [figure_dir, 'microstructure.png']);
end

function build_graph_scatters(GM,data,modality_name, figure_dir)
% Build microstructure 2D scatter plots.
%
%   BUILD_GRAPH_SCATTERS(GM,data,modality_name,figure_dir) builds scatter
%   plots of eccentricity versus vectors in data. Scatter plots are stored
%   in figure_dir/{side}_{modality}.png where side is left/right and
%   modality is curvature, cortical thickness or t1wt2w. 

hemi_name = ["lh", "rh"];
import temporal_gradients.support.eccentricity
import temporal_gradients.support.metric_scatter

for hemi = 1:2
    for modality = 1:size(data.(hemi_name(hemi)),2)
        if modality == 1
            ylim = [-0.3 0.4];
            ylab = 'Curvature';
        elseif modality == 2
            ylim = [1.5 4.5];
            ylab = 'Cortical Thickness';
        elseif modality == 3
            ylim = [0.7 2.5];
            ylab = 'T1w/T2w';
        end
        name = figure_dir + ...
            "scatter_" + ...
            hemi_name(hemi) + "_" + modality_name(modality) +".png";
        metric_scatter(eccentricity(GM.aligned{hemi}(:,1:3)), ...
            data.(hemi_name{hemi})(:,modality), ...
            [0,5],ylim,'Eccentricity',ylab,true,replace(name,'T1w/T2w.png','T1wT2w.png'));
    end
end
end


function build_graph_scatters_3d(GM,data,clim,cmap,modality_name, figure_dir)
% Build microstructure 3D scatter plots.
%
%   BUILD_GRAPH_SCATTERS_3D(GM,data,clim,cmap,modality_name,figure_dir)
%   builds scatter plots of vectors in data in the manifold space described
%   by the first three gradients in the GradientMaps object GM. Variables
%   clim and cmap provide the color limits and colormap, respectively.
%   Scatter plots are stored in figure_dir/{side}_{modality}.png where side
%   is left/right and modality is derived cell/string array modality_name.

hemi_name = ["lh", "rh"];
import temporal_gradients.support.eccentricity
import temporal_gradients.support.metric_scatter_3d

for hemi = 1:2
    for modality = 1:size(data.(hemi_name(hemi)),2)
        name = figure_dir + ...
            "scatter3d_" + ...
            hemi_name(hemi) + "_" + modality_name(modality) +".png";
        metric_scatter_3d(GM.aligned{hemi}, data.(hemi_name{hemi})(:,modality), ...
            clim(modality,:),replace(name,'T1w/T2w.png','T1wT2w.png'), cmap{modality});
    end
end
end


function build_graph_bigbrain(GM,data,parcellation, figure_dir)
% Build bigbrain scatter scatter plots.
%
%   BUILD_GRAPH_BIGBRAIN(GM,data,parcellaton,figure_dir) builds scatter
%   plots of eccentricity derived from the GradientMaps object GM versus
%   the vector data. Eccentricity is parcellated using the parcellation
%   vector. The figure is stored as
%   figure_dir/scatter_{side}_mpc_bigbrain.png.

hemi_name = ["lh", "rh"];
import temporal_gradients.support.eccentricity
import temporal_gradients.support.metric_scatter

for hemi = 1:2
    if hemi == 1 
        p = parcellation(1:end/2);
        d = labelmean(data(1:end/2),p);
    else
        p = parcellation(end/2+1:end);
        d = labelmean(data(end/2+1:end),p);
    end
    d(isnan(d)) = [];
    downsample = labelmean(eccentricity(GM.aligned{hemi}(:,1:3)),p);
    downsample(isnan(downsample)) = []; 
    name = figure_dir + ...
        "scatter_" + ...
        hemi_name(hemi) + "_mpc_bigbrain.png";
    metric_scatter(downsample, ...
                   d, ...
                   [0,4],[-.15 .2],'Eccentricity','MPC',true,name);
end
end

function build_surfaces_tl(data,surfaces,clim,cmap,labels,mask,figure_name)
% Builds cortical surfaces
%
% 	BUILD_SURFACES(data,surfaces,clim,cmap,labels,mask,figure_name)
% 	visualizes columns of matrix data on the cortical surfaces with color
% 	limits and maps clim and cmap. A label is plotted for each data vector
% 	from the labels in string/cell array labels. A mask is provided to
% 	denote the vertices included in the temporal lobe. The figure is stored
% 	as figure_name.

fake_parcellation = zeros(20000,1);
fake_parcellation(mask) = 1:3428;
obj = plot_hemispheres(data,surfaces,'parcellation',fake_parcellation','labeltext',labels);
obj.colorlimits(clim);
obj.colormaps(cmap);
export_fig(figure_name,'-m2','-png');
end

function cmap = interp_rgb(val,sz)
% Constructs a colormap by linear interpolation.
%
% cmap = INTERP_RGB(val,sz) returns a sz-by-3 matrix of values interpolated
% between the rows of N-by-3 matrix val. 

x = linspace(0,1,size(val,1)); 
cmap = zeros(sz,3);
for ii = 1:3
    cmap(end:-1:1,ii) = interp1(x,val(:,ii),linspace(0,1,sz))/255;
end
end

function cmap = greensequential(sz)
% Constructs a green sequential colormap.
%
%   cmap = GREENSEQUENTIAL(sz) constructs a green sequential colormap of
%   size sz-by-3.

if nargin == 0
    sz = 64;
end
cmap = interp_rgb([247,252,245
                   229,245,224
                   199,233,192
                   161,217,155
                   116,196,118
                   65,171,93
                   35,139,69
                   0,109,44
                   0,68,27],sz);
end