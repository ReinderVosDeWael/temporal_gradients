function download_data()
% Downloads data required for figure generation
%
%   DOWNLOAD_DATA() downloads data required for building the figures for
%   the Vos de Wael et al., 2020, bioRxiv. The data is stored inside the
%   data directory within this package. Please note, the downloader looks
%   for the package directory by name, i.e. +temporal_gradients, and so
%   will not work if this directory is renamed. 
%
%   Written by Reinder Vos de Wael, MICA lab, Nov 2020
%   For further details see our <a
%   href="https://github.com/MICA-MNI/micaopen/tree/master/temporal_gradients">Github</a> page.  

fs = string(filesep);
package_dir = regexp(mfilename('fullpath'),'.*\+temporal_gradients','match','once');
websave(package_dir + fs + 'data' + fs + 'figure_data.mat', 'https://box.bic.mni.mcgill.ca/s/M0WDuW61Q3eFL3t/download');
end