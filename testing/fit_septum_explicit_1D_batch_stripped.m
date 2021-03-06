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

function ud = fit_septum_explicit_1D_batch_stripped(fullstack, tracks, param, xy0, olaystack, tracknums)

if nargin < 6
    tracknums = 1:param.ntracks;
end
if nargin < 4
    xy0 = [];
    olaystack = [];
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
        fname_end = param.im_file(end-6:end);
        fname = [param.path '/' param.analysis_date '_' fname_base '_pos' fovnum '_const_raw.fig'];
    end
    
    gr = figure('FileName', fname, 'Position', [100 100 500 400]);
%     gr = figure('FileName', [param.path '/' param.analysis_date]);
    hr = subplot(211);
    box on
    hold on
    title(['Raw constriction traces for ' strrep(param.im_file,'_','\_')])
    ylabel('Diameter (nm)')
    
    hh = subplot(212);
    hold on
    box on
    ylabel('Thickness (nm)')
%     subplot(3,1,1);
%     subplot(3,1,2);
%     ylabel('Axial width (nm)')
%     subplot(3,1,3);
%     ylabel('Intensity')
    xlabel('Time (min)')
end

for ii = tracknums
    %% Load image stack and pull out individual septa
    
    display(num2str(ii))
    if ~isfield(tracks.Experiment.Lineage(ii).Trajectory(1).LIFESPAN,'f')
        disp('trajectory removed in microbeJ.')
        continue
    elseif length(tracks.Experiment.Lineage(ii).Trajectory)>1
        disp('microbeJ recorded division event. removing.')
        continue
    else
        [stack, centx, centy, imframes, olaystack_cropped, xy0_subim] = separate_microbej_tracks(fullstack, tracks, ii, param.s_box, param.n_frames_before, olaystack);
    end
    
    %% Obtain line profiles for all frames
    
    [improfs_rad, improfs_ax, fit_rad, fit_ax] = fit_septum_explicit_1D_general(stack, param.plot_im, param, xy0_subim, olaystack_cropped);

    avail_frames_before = imframes(1); % how many frames before can we use?
    frames_before = min([param.n_frames_before avail_frames_before]);
    start_frame_new = max([imframes(1) - param.n_frames_before, 1]); % add a few extra frames before
    frames = ([(start_frame_new:start_frame_new+frames_before-1) imframes])'; % removed + 1 200124
    
    frames = double(frames);
    
    ud.imDat(ii).num = ii;
    ud.imDat(ii).line_profiles = improfs_rad;
    ud.imDat(ii).frames = (imframes + 1)';
    ud.imDat(ii).centx = centx';
    ud.imDat(ii).centy = centy';
    
    % add in excluded frames (that were out of focus or whatever)
    if ~isempty(param.exclude_frames) && any(frames(1)<param.exclude_frames) && any(frames(end)>param.exclude_frames)
        xx = param.exclude_frames;
        [exind, ~] = find(frames==xx);
        for kk = 1:length(exind)
            frames = [frames(1:exind(kk)-1); xx(kk); frames(exind(kk):end)+1];
            %         centx = [centx(1:exind-1) NaN centx(exind:end)];
            %         centy = [centy(1:exind-1) NaN centy(exind:end)];
            improfs_rad = [improfs_rad(1:exind(kk)-1,:); ones(1,size(improfs_rad,2))*NaN; improfs_rad(exind(kk):end,:)];
            
            improfs_ax = [improfs_ax(1:exind(kk)-1,:); ones(1,size(improfs_ax,2))*NaN; improfs_ax(exind(kk):end,:)];

            fit_rad = [fit_rad(1:exind(kk)-1,:); ones(1,size(fit_rad,2))*NaN; fit_rad(exind(kk):end,:)];
            fit_ax = [fit_ax(1:exind(kk)-1,:); ones(1,size(fit_ax,2))*NaN; fit_ax(exind(kk):end,:)];
        end
    end
    
    linecents = fit_rad(:,1); % [pix]
    diams = abs(fit_rad(:,2)) *2*param.pixSz; % [nm]
    fiterrs = fit_rad(:,3); % [intensity^2] residual sum of squares from fit to line profile
    
    FWHM_ax = fit_ax(:,1) *param.pixSz; % [nm]
    se_ax = fit_ax(:,2); % should probably be *2*param.pixSz
    
%     FWHM_ax = fit_ax(:,2) *2*param.pixSz; % [nm]
%     se_ax = fit_ax(:,3); % [intensity^2]
    
    time = frames * param.interval; % [min]
    
    %% Fit line profiles to generalized Gaussian function to get FWHM
    
    intensity = [];
    for jj = 1:size(improfs_rad,1)
        linep = improfs_rad(jj,~isnan(improfs_rad(jj,:)));
        if isempty(linep) || length(linep)<5
            intensity(jj) = NaN;
            intensityMedian(jj) = NaN;
        else
            intensity(jj) = sum(linep);
            intensityMedian(jj) = median(linep);
        end
    end

    %start with vars
    %time, intensity,se_ax,fiterrs (radial SE), diams (radial width),
    %FWHM_ax (axial width)
    %ud.rawDat(ii).cuttime = tCrop_rad;
    %make a struct array ud.datUnfilt
    ud.datUnfilt(ii).time = time;
    ud.datUnfilt(ii).intensity = intensity;
    ud.datUnfilt(ii).intensityMedian = intensityMedian;
    ud.datUnfilt(ii).se_ax = se_ax;
    ud.datUnfilt(ii).fiterrs = fiterrs;
    ud.datUnfilt(ii).linecents = linecents;
    ud.datUnfilt(ii).diam = diams;
    ud.datUnfilt(ii).FWHM_ax = FWHM_ax;
    
    %% Filter data

    MAXDIAM = 1400;
    MINDIAM = 50;
%     FWHM_ax = FWHM_ax.*param.pixSz;
    t_ax = time;
    t_rad = time;
    int_ax = intensity;
    int_rad = intensity;
    %indepently each axis filters out too wide, too narrow
    badIdx_ax = FWHM_ax<MINDIAM | isnan(FWHM_ax);
    badIdx_rad = diams<MINDIAM | diams>MAXDIAM | isnan(diams);
%     FWHM_ax(badIdx_ax) = [];
%     t_ax(badIdx_ax) = [];
%     se_ax(badIdx_ax) = [];
%     int_ax(badIdx_ax) = [];
    FWHM_ax(badIdx_rad) = [];
    t_ax(badIdx_rad) = [];
    se_ax(badIdx_rad) = [];
    int_ax(badIdx_rad) = [];
    diams(badIdx_rad) = [];
    linecents(badIdx_rad) = [];
    t_rad(badIdx_rad) = [];
    int_rad(badIdx_rad) = [];
    fiterrs(badIdx_rad) = [];

    %independently each axis filters out too big fitting error
%     badIdx3_ax = se_ax(4,:) > 2;
%     badIdx3_ax = se_ax > 0.5;
    badIdx3_ax = se_ax > 3e3; % SSR
    badIdx3_ax = se_ax < 0.9; % R^2
%     FWHM_ax(badIdx3_ax) = [];
%     t_ax(badIdx3_ax) = [];
%     int_ax(badIdx3_ax) = [];
%     badIdx3_rad = fiterrs > 1.1;
    badIdx3_rad = fiterrs > 2e4; % SSR
    badIdx3_rad = fiterrs < 0.9; % R^2
    diams(badIdx3_rad) = [];
    linecents(badIdx3_rad) = [];
    t_rad(badIdx3_rad) = [];
    int_rad(badIdx3_rad) = [];
    FWHM_ax(badIdx3_rad) = [];
    t_ax(badIdx3_rad) = [];
    int_ax(badIdx3_rad) = [];
    
    % filter based on fitted center of line profile - if very off from
    % center (relative to fitted diameter) it's likely a poor fit even with
    % low fit error.
    badIdx4_rad = linecents<diams/param.pixSz*2/3 | size(improfs_rad,2)-linecents<diams/param.pixSz*2/3;
    diams(badIdx4_rad) = [];
    linecents(badIdx4_rad) = [];
    t_rad(badIdx4_rad) = [];
    int_rad(badIdx4_rad) = [];
    
    ud.rawDat(ii).num = ii;
    ud.rawDat(ii).width = diams;
    ud.rawDat(ii).linecents = linecents;
    ud.rawDat(ii).time = t_rad;
    ud.rawDat(ii).intensity = int_rad; % added 200120
    ud.rawDat(ii).width_ax = FWHM_ax;
    ud.rawDat(ii).time_ax = t_ax;
    
    % need a way to delete all the crap after constriction finishes
    % Constriction end is just after intensity peak
    if param.intensity_cut
        [iMax_ax, idxMax_ax] = max(int_ax);
        [iMax_rad, idxMax_ex] = max(int_rad);

        iMin_ax = min(int_ax(t_ax>mean(t_ax)));
        iMin_rad = min(int_rad(t_rad>mean(t_rad)));

        iConsEnd_ax = (iMax_ax+iMin_ax)/2;
        tCrop_ax = t_ax(idxMax_ax:end);
        iCrop_ax = int_ax(idxMax_ax:end);
%         iConsEnd_rad = (iMax_rad+iMin_rad)/2;
        iConsEnd_rad = (iMax_rad-iMin_rad)*0.8 + iMin_rad; % 80% of max
        tCrop_rad = t_rad(idxMax_ex:end);
        iCrop_rad = int_rad(idxMax_ex:end);

        iDiff_ax = iCrop_ax - iConsEnd_ax;
        tConsEnd_ax = tCrop_ax(find(iDiff_ax<0,1)-1);
        iDiff_rad = iCrop_rad - iConsEnd_rad;
        tConsEnd_rad = tCrop_rad(find(iDiff_rad<0,1)-1);

        if isempty(tConsEnd_ax)
            tidx_ax = length(t_ax);
        else
            tidx_ax = find(t_ax==tConsEnd_ax);
        end
        if isempty(tConsEnd_rad)
            tidx_rad = length(t_rad);
        else
            tidx_rad = find(t_rad==tConsEnd_rad);
        end

        dCrop = diams(1:tidx_rad);
        cCrop = linecents(1:tidx_rad);
        iCrop = int_rad(1:tidx_rad);
        tCrop_ax = t_ax(1:tidx_ax);
        fCrop_ax = FWHM_ax(1:tidx_ax);
        tCrop_rad = t_rad(1:tidx_rad);
    else
        dCrop = diams;
        cCrop = linecents;
        iCrop = int_rad;
        tCrop_ax = t_ax;
        fCrop_ax = FWHM_ax;
        tCrop_rad = t_rad;
    end
    
    ud.rawDat(ii).cuttime = tCrop_rad;
    ud.rawDat(ii).cutdiams = dCrop;
    ud.rawDat(ii).cutcents = cCrop;
    ud.rawDat(ii).cutint = iCrop;
    ud.rawDat(ii).cuttime_ax = tCrop_ax;
    ud.rawDat(ii).cutdiams_ax = fCrop_ax;

    if param.plot_raw
%         plot(hr, tCrop_rad, dCrop)
        
%         figure(hh)
%         subplot(2,1,1)
        plot(hr, tCrop_rad, dCrop)
%         subplot(2,1,2)
        plot(hh, tCrop_ax,fCrop_ax)
%         subplot(3,1,3);
%         plot(ud.rawDat(ii).time,ud.rawDat(ii).intensity)

    end

end

ud.param = param;

if param.plot_raw
    set(gr, 'UserData', ud)
    
    if param.save_raw
        savefig(gr, fname)
    end
end
