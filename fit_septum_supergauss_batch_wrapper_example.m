% Author: Kevin Whitley
% Date created: 190729

% This script asks the user for a directory containing a time-lapse image
% file (.tif) and a microbeJ track output (.mat). It contains all the user
% input parameters that will later be saved in the figure file. It then
% calls fit_septum_supergauss_batch_stripped, which will in turn call
% separate_microbej_tracks, fit_septum_supergauss, and
% fit_constriction_data.

%% User-set parameters

param.plot_im = 0; % plot septum images with fitted diameter
param.plot_raw = 1; % plot raw constriction vs. time data
param.plot_filt = 0; % filtered data (crap removed)

% params for separate_microbej_tracks
param.s_box = 24; % size of sub-images
param.n_frames_before = 0; % number of frames before first one to include

% params for fit_septum_supergauss
param.psfFWHM = 250; % [nm]
param.pixSz = 65; % [nm/pix] Theia

param.interval = 1; % [min/frame]

param.t_cpd = 0; % [min] time at which compound flushed in

param.exclude_frames = []; % frames where focus was lost, etc. need to add NaNs to everything there.

%% Files

today = datestr(now, 'yymmdd');

path = pwd;
dirIm = dir([path '\191128_3_MMStack_Pos2.ome_denoise_reg_cut_bgsub.tif']);
dirTr = dir([path '\200108_191128_3_pos2_microbej.mat']);

param.path = path;
param.analysis_date = today;
param.im_file = dirIm.name;
param.tracks_file = dirTr.name;

fullstack = imreadstack([path '\' dirIm.name]);
fullstack = im2double(fullstack);

tracks = load([path '\' dirTr.name]);
param.ntracks = length(tracks.Experiment.Lineage);

%% Run function

[ud, t_con] = fit_septum_supergauss_batch_stripped(fullstack, tracks, param);



