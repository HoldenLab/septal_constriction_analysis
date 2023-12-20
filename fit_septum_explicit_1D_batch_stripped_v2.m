% Author: Kevin Whitley
% Date created: 191007

% This function calls several functions to fit bacterial septa to an
% explicit model and analyze the resulting trace of diameter vs. time:

% 1. separate_microbej_tracks cuts one tracked septum out of the original
% image with a small bounding window.

% 2. fit_septum_explicit takes this time-lapse of one septum and fits each
% frame to an explicit septum-image model.

% 3. fit_constriction_data fits the septum diameter vs. time to a model of
% the user's choice.

% Output variable ud (for 'user data') is a structure variable with two
% fields:
%   - imDat: contains frames, centroid positions, and line profiles for
%   each individual fitted image
%   - unfiltDat: contains unfiltered results of fitting each line profile
%   to its appropriate model (diameter, thickness, fitting errors, etc.)

% Output variable ud is added to the UserData of the output figure. So, all
% resulting data can be saved and accessed through the figure itself - pull
% up figure, then use [variable name] = get(gcf,'UserData').

% v2 (201202):
%   - Removed all filtering and field 'rawDat' from output variable ud. All
%   filtering now done in plot_constriction_results2.
%   - Removed option to exclude out of focus frames. No longer necessary,
%   since any very poor fits will be filtered out later anyway.

function ud = fit_septum_explicit_1D_batch_stripped_v2(fullstack, tracks, param, olaystack, tracknums)

if nargin < 5
    tracknums = 1:param.ntracks;
end

if param.plot_raw
    
    % figure out file name format (different with different software)
    isMM = strfind(param.im_file, '.ome');
    if ~isempty(isMM) % micro-manager format
        fname_base = extractBefore(param.im_file, '.ome');
        fname_base = erase(fname_base, 'MMStack_');
        fname = [param.path '/' param.analysis_date '_' fname_base '_const_raw.fig'];
    else % NS elements format
        fname_base = extractBefore(param.im_file, '.nd2 -');
        fovstr = extractAfter(param.im_file, '.nd2 (series ');
        fovnum = fovstr(1:2);
        fname = [param.path '/' param.analysis_date '_' fname_base '_pos' fovnum '_const_raw.fig'];
    end
    
    % set up figure for diameters/thicknesses
    h_fig = figure('FileName', fname, 'Position', [100 100 500 400]);
    
    h_diam = subplot(311, 'box', 'on');
    hold on
    title(['Raw constriction traces for ' strrep(param.im_file,'_','\_')])
    ylabel('Diameter (nm)')
    
    h_thick = subplot(312, 'box', 'on');
    hold on
    ylabel('Thickness (nm)')
    
    h_int = subplot(313, 'box', 'on');
    hold on
    ylabel('Intensity')
    xlabel('Time (min)')
    
end

for ii = tracknums
    %% Load image stack and pull out individual septa
    
    display(num2str(ii))
    
    if ~isfield(tracks.Experiment.Lineage(ii).Trajectory(1).LIFESPAN,'f')
        disp('trajectory removed in microbeJ.')
        continue
    elseif length(tracks.Experiment.Lineage(ii).Trajectory) > 1
        disp('microbeJ recorded division event. removing.')
        continue
    else
        [stack, centx, centy, imframes, olaystack_cropped, xy0_subim] = separate_microbej_tracks_v2(fullstack, tracks, ii, param, olaystack);
    end
    
    %% Obtain and fit line profiles for all frames
    
    [improfs_rad, improfs_ax, fit_rad, fit_ax] = fit_septum_explicit_1D_general_v2(stack, param, xy0_subim, olaystack_cropped);

    % add extra frames before beginning of tracks to match line profiles (if requested by user)
    avail_frames_before = imframes(1); % how many frames before can we use?
    frames_before = min([param.n_frames_before avail_frames_before]);
    start_frame_new = max([imframes(1) - param.n_frames_before, 1]); % add a few extra frames before
    frames = ([(start_frame_new:start_frame_new+frames_before-1) imframes])'; % removed + 1 200124
    frames = double(frames); % necessary?
    
    ud.imDat(ii).num = ii;
    ud.imDat(ii).frames = (imframes + 1);
    ud.imDat(ii).centx = centx;
    ud.imDat(ii).centy = centy;
    ud.imDat(ii).line_profiles = improfs_rad;
    ud.imDat(ii).ax_line_profiles = improfs_ax;

    %% Grab unfiltered data
    
    time = frames * param.interval; % [min]
    
    linecents_rad = fit_rad(:,1); % [pix] centers of (radial) line profiles
    diams = abs(fit_rad(:,2)) *2*param.pixSz; % [nm]
    fiterrs_rad = fit_rad(:,3); % [intensity^2] residual sum of squares from fit to line profile
    
    FWHM_ax = fit_ax(:,1) *param.pixSz; % [nm]
    se_ax = fit_ax(:,2); % R^2 from fit to super-gaussian
    std_ax = fit_ax(:,3);

    intensity = nansum(improfs_rad, 2);
    intensityMedian = median(improfs_rad, 2, 'omitnan');
    
    int_ex = sum(~isnan(improfs_rad),2) <= 5; % profiles with fewer than 5 real values
    intensity(int_ex) = NaN; % exclude if too few real values
    intensityMedian(int_ex) = NaN;
    
    int_ax = nansum(improfs_ax, 2);

    ud.datUnfilt(ii).num = ii;
    ud.datUnfilt(ii).time = time;
    ud.datUnfilt(ii).linecents = linecents_rad;
    ud.datUnfilt(ii).diam = diams;
    ud.datUnfilt(ii).fiterrs = fiterrs_rad;
    ud.datUnfilt(ii).FWHM_ax = FWHM_ax;
    ud.datUnfilt(ii).se_ax = se_ax;
    ud.datUnfilt(ii).int_ax = int_ax;
    ud.datUnfilt(ii).std_ax = std_ax;
    ud.datUnfilt(ii).intensity = intensity;
    ud.datUnfilt(ii).intensityMedian = intensityMedian;

    if param.plot_raw

        plot(h_diam, time, diams)
        plot(h_thick, time, FWHM_ax)
        plot(h_int, time, intensity)

    end

end

ud.param = param;

if param.plot_raw
    set(h_fig, 'UserData', ud) % add user data (output variable) to figure
    
    if param.save_raw
        savefig(h_fig, fname)
    end
end
