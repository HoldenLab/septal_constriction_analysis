% Author: Kevin Whitley
% Date created: 190729

% This script asks the user for a directory containing a time-lapse image
% file (.tif) and a microbeJ track output (.mat). It contains all the user
% input parameters that will later be saved in the figure file. It then
% calls fit_septum_explicit_batch_stripped, which will in turn call
% separate_microbej_tracks, fit_septum_explicit, and fit_constriction_data.

%% User-set parameters

param.ZGFP = 1; % are these images ZGFP strain? changes some image processing. 
%CURRENTLY ZGFP uses the olay image
param.plot_im = 0; % plot septum images with fitted diameter
param.plot_gauss = 0;
param.plot_explicit = 0;
param.plot_raw = 0; % plot raw constriction vs. time data
param.save_raw = 0;
param.intensity_cut = 0;

% params for separate_microbej_tracks
param.s_box = 30; % size of sub-images
param.n_frames_before = 0; % number of frames before first one to include

% params for fit_septum_explicit
param.psfFWHM = 250; % [nm]
param.pixSz = 65; % [nm/pix] Theia

param.interval = 1; % [min/frame]

param.t_cpd = 0; % [min] time at which compound flushed in

param.exclude_frames = []; % frames where focus was lost, etc. need to add NaNs to everything there.

%% Files

today = datestr(now, 'yymmdd');

path = pwd;
dirIm = dir([path '\200104_1mW-TIRF-30_CDM_1_MMStack_Default.ome_denoise_reg_32bit.tif']);
if param.ZGFP
    dirImOlay = dir([path '\200104_1mW-TIRF-30_CDM_1_MMStack_Default.ome_denoise_reg_32bit_Simple Segmentation_.tif']);
end
dirTr = dir([path '\Segmented Image Trial analysis.mat']);

param.path = path;
param.analysis_date = today;
param.im_file = dirIm.name;
param.tracks_file = dirTr.name;

fullstack = imreadstack([path '\' dirIm.name]);
fullstack = im2double(fullstack);
if param.ZGFP
    olaystack = imreadstack([path '\' dirImOlay.name]);
    olaystack = im2double(olaystack);
end

tracks = load([path '\' dirTr.name]);
param.ntracks = length(tracks.Experiment.Lineage);
%% Run function

if param.ZGFP
    %for ii =1:199
    %    close all
    %    ud = fit_septum_explicit_1D_batch_stripped(fullstack, tracks, param,ii,olaystack);
    %end
    
    %ud = fit_septum_explicit_1D_batch_stripped(fullstack, tracks, param,[20 47],olaystack);
    %ud = fit_septum_explicit_1D_batch_stripped(fullstack, tracks, param,[39],olaystack);
    ud = fit_septum_explicit_1D_batch_stripped(fullstack, tracks, param,[],olaystack);
else
    ud = fit_septum_explicit_1D_batch_stripped(fullstack, tracks, param);
end

