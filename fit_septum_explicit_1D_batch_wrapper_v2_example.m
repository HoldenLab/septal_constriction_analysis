% Author: Kevin Whitley
% Date created: 190729

% This script asks the user for a directory containing a time-lapse image
% file (.tif) and a microbeJ track output (.mat). It contains all the user
% input parameters that will later be saved in the figure file. It then
% calls fit_septum_explicit_batch_stripped, which will in turn call
% separate_microbej_tracks, fit_septum_explicit, and fit_constriction_data.

%% User-set parameters

param.ZGFP = 1; % are these images ZGFP strain? changes one small part of image processing if not doing segmentation
param.segment = 1; % use image stack with segmented data from ilastik

param.plot_im_fits = 0; % plot septum images for each frame along with line profiles of diameter and thickness with model fits. for debugging.
param.plot_raw = 1; % plot septal diameter and thickness vs. time
param.save_raw = 1; % automatically save fig file with all results contained in 'UserData'

% params for separate_microbej_tracks
param.s_box = 30; % size of sub-images
param.n_frames_before = 0; % number of frames before first one to include

% microscope parameters
param.psfFWHM = 250; % [nm]
param.pixSz = 65; % [nm/pix] Theia

param.interval = 1; % [min/frame]

param.t_cpd = 0; % [min] time at which compound flushed in

%% Load files

path = pwd;

% load image stack
dirIm = dir([path '\*bgsub.tif']);
fullstack = imreadstack([path '\' dirIm.name]);
fullstack = im2double(fullstack);

% load tracking file (microbej matlab output)
dirTr = dir([path '\*microbej.mat']);
tracks = load([path '\' dirTr.name]);

% load image stack for segmented data (from ilastik) if doing ZGFP strain
if param.segment
    dirImOlay = dir([path '\segmented\*binary.tif']);
    olaystack = imreadstack([path '\segmented\' dirImOlay.name]);
    olaystack = im2double(olaystack);
else
    olaystack = [];
end

param.path = path;
param.analysis_date = datestr(now, 'yymmdd'); % today's date
param.im_file = dirIm.name;
param.tracks_file = dirTr.name;
param.ntracks = length(tracks.Experiment.Lineage);

%% Run function

ud = fit_septum_explicit_1D_batch_stripped_v2(fullstack, tracks, param, olaystack);

