% Author: Kevin Whitley
% Date created: 190729

% This script asks the user for a directory containing a time-lapse image
% file (.tif) and a microbeJ track output (.mat). It contains all the user
% input parameters that will later be saved in the figure file. It then
% calls fit_septum_explicit_batch_stripped, which will in turn call
% separate_microbej_tracks, fit_septum_explicit, and fit_constriction_data.

%% User-set parameters

param.plot_im = 0; % plot septum images with fitted diameter
param.plot_raw = 0; % plot raw constriction vs. time data
param.plot_filt = 1; % filtered data (crap removed)
plot_t_con = 0; % histogram of constriction times

% params for separate_microbej_tracks
param.s_box = 24; % size of sub-images
param.n_frames_before = 5; % number of frames before first one to include

% params for fit_septum_explicit
param.psfFWHM = 250; % [nm]
param.pixSz = 65; % [nm/pix] Theia
% param.pixSz = 106; % [nm/pix] NSIM
param.gm_model = 0; % use mixed Gaussian model (else - tilted circle model)
param.gridsp = 0.5;
param.fval_thresh = 0.006;
param.cc_thresh = 0.77;
param.ssim_thresh = 0.94;

orientation_cut = 0; % seems to not work so well. leave off for now.
param.interval = 1; % [min/frame]
param.fit_model = 1;

% params for fit_constriction_data
param.model = 'parabolic';
param.chi_thresh = 2e4;

param.t_cpd = 38; % [min] time at which compound flushed in

param.exclude_frames = []; % frames where focus was lost, etc. need to add NaNs to everything there.

%% Files

today = datestr(now, 'yymmdd');

path = uigetdir('D:\Data_Theia\');
dirIm = dir([path '\*8bit.tif']);
dirTr = dir([path '\*microbej.mat']);

param.path = path;
param.analysis_date = today;
param.im_file = dirIm.name;
param.tracks_file = dirTr.name;

fullstack = imreadstack([path '\' dirIm.name]);
fullstack = im2double(fullstack);

tracks = load([path '\' dirTr.name]);
param.ntracks = length(tracks.Experiment.Lineage);

%% Run function

[ud, t_con] = fit_septum_explicit_batch_stripped(fullstack, tracks, param);

if plot_t_con
    figure('FileName', [path '/' today '_t_con.fig'])
    t_con(t_con<1) = [];
    hist(t_con)
    xlabel('t_{constriction} (min)')
    ylabel('Counts')
end
