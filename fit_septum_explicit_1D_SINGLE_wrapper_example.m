% Author: Kevin Whitley
% Date created: 200128

% This script asks the user for a directory containing a time-lapse image
% file (.tif) and a microbeJ track output (.mat). It contains all the user
% input parameters that will later be saved in the figure file. It then
% calls fit_septum_explicit_batch_stripped, which will in turn call
% separate_microbej_tracks, fit_septum_explicit, and fit_constriction_data.

%% User-set parameters

param.ZGFP = 1;
param.segment = 1; % use image stack with segmented data from ilastik

param.plot_im_fits = 0; % plot septum images for each frame along with line profiles of diameter and thickness with model fits. for debugging.
param.plot_raw = 0; % plot raw constriction vs. time data

% params for separate_microbej_tracks
param.s_box = 30; % size of sub-images
param.n_frames_before = 0; % number of frames before first one to include

% params for fit_septum_explicit
param.psfFWHM = 250; % [nm]
param.pixSz = 65; % [nm/pix] Theia

param.interval = 1; % [min/frame]

param.t_cpd = 33; % [min] time at which compound flushed in

%% Files

today = datestr(now, 'yymmdd');

path = [pwd '\Indiv_rings'];
dirIm = dir([path '\*.tif']);
dirTr = [];

if param.ZGFP
    olaypath = [pwd '\segmented\Indiv_rings'];
    dirImOlay = dir([olaypath '\*.tif']);
end

param.path = path;
param.analysis_date = today;

fits = []; results = []; fitex = [];
for ii = 1:length(dirIm)
    
    param.im_file = dirIm(ii).name;
    param.tracks_file = [];
    
    imstack = imreadstack([path '\' dirIm(ii).name]);
    imstack = im2double(imstack);
    
    if param.segment
        olaystack = imreadstack([olaypath '\' dirImOlay(ii).name]);
        olaystack = im2double(olaystack);
    end
    
    param.ntracks = 0;
    
    param.date = str2double(param.im_file(1:6));
    param.filenum = str2double(param.im_file(8));
    posind = strfind(param.im_file, 'Pos');
    param.pos = str2double(param.im_file(posind+3));
    ringind = strfind(param.im_file, 'ring');
    im_file_trunc = param.im_file(1:end-4); % removes .tif from end
    param.ringnum = str2double(im_file_trunc(ringind+4:end));
    
    %% Run function

    [improf, orthprof, fitex, fitax] = fit_septum_explicit_1D_general_v2(imstack, param, [14 14], olaystack);

    results(ii,:) = [param.date param.filenum param.pos param.ringnum fitex fitax];
    
end

results = sortrows(results, 4);

