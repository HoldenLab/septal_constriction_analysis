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

function ud = fit_septum_explicit_1D_batch_stripped(fullstack, tracks, param, tracknums)

if nargin < 4
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
        fname_end = param.im_file(end-6:end);
        fname = [param.path '/' param.analysis_date '_' fname_base '_pos' fname_end '_const_raw.fig'];
%         fname = [param.path '/' param.analysis_date '_.fig'];
    end
    
    gr = figure('FileName', fname, 'Position', [100 100 500 400]);
%     gr = figure('FileName', [param.path '/' param.analysis_date]);
    hh = gca;
    hold on
    box on
    title(['Raw constriction traces for ' strrep(param.im_file,'_','\_')])
    xlabel(hh,'Time (min)')
    ylabel(hh,'Septum width (nm)')
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
        [stack, centx, centy, imframes] = separate_microbej_tracks(fullstack, tracks, ii, param.s_box, param.n_frames_before);
    end
    
    %% Obtain line profiles for all frames
    
%     [improfs_rad, improfs_ax, fit_rad, fit_ax] = fit_septum_explicit_1D(stack, param.plot_im, param);
    [improfs_rad, improfs_ax, fit_rad, fit_ax] = fit_septum_explicit_1D_imrot(stack, param.plot_im, param);

    avail_frames_before = imframes(1); % how many frames before can we use?
    frames_before = min([param.n_frames_before avail_frames_before]);
    start_frame_new = max([imframes(1) - param.n_frames_before, 1]); % add a few extra frames before
    frames = ([(start_frame_new:start_frame_new+frames_before-1) imframes])'; % removed + 1 200124
%     frames = ([(start_frame_new:start_frame_new+frames_before-1) imframes] + 1)';
    
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
    
    diams = abs(fit_rad(:,2)) *2*param.pixSz;
    fiterrs = fit_rad(:,3); % should probably be *2*param.pixSz
    
    FWHM_ax = fit_ax(:,1) *2*param.pixSz;
    se_ax = fit_ax(:,2); % should probably be *2*param.pixSz
    
    time = frames * param.interval; % [min]
    
    %% Fit line profiles to generalized Gaussian function to get FWHM
    
    intensity = [];
    for jj = 1:size(improfs_rad,1)
        linep = improfs_rad(jj,~isnan(improfs_rad(jj,:)));
        if isempty(linep) || length(linep)<5
            intensity(jj) = NaN;
        else
%             [~, intensity(jj)] = fitProfile_public(linep,[],1);
            intensity(jj) = sum(linep);
        end
    end
%     for kk = 1:size(improfs_ax,1)
%         linep_ax = improfs_ax(kk,~isnan(improfs_ax(kk,:)));
%         if isempty(linep_ax) || length(linep_ax)<5
%             FWHM_ax(kk) = NaN;
%             se_ax(:,kk) = ones(1,5)*NaN;
%         else
%             [FWHM_ax(kk), ~, ~, ~, se_ax(:,kk)] = fitProfile_public(linep_ax,[],1); % axial FWHM
%         end
%     end

    %% Cut results

    MAXDIAM = 1400;
    MINDIAM = 50;
%     FWHM_ax = FWHM_ax.*param.pixSz;
    t_ax = time;
    t_rad = time;
    int_ax = intensity;
    int_rad = intensity;
    badIdx_ax = FWHM_ax<MINDIAM | isnan(FWHM_ax);
    badIdx_rad = diams<MINDIAM | diams>MAXDIAM | isnan(diams);
    FWHM_ax(badIdx_ax) = [];
    t_ax(badIdx_ax) = [];
    se_ax(badIdx_ax) = [];
    int_ax(badIdx_ax) = [];
    diams(badIdx_rad) = [];
    t_rad(badIdx_rad) = [];
    int_rad(badIdx_rad) = [];
    fiterrs(badIdx_rad) = [];

%     badIdx3_ax = se_ax(4,:) > 2;
    badIdx3_ax = se_ax > 0.5;
    FWHM_ax(badIdx3_ax) = [];
    t_ax(badIdx3_ax) = [];
    int_ax(badIdx3_ax) = [];
    badIdx3_rad = fiterrs > 1.1;
%     badIdx3_rad = fiterrs > 0.5;
    diams(badIdx3_rad) = [];
    t_rad(badIdx3_rad) = [];
    int_rad(badIdx3_rad) = [];
    
    ud.rawDat(ii).num = ii;
    ud.rawDat(ii).width = diams;
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
        iConsEnd_rad = (iMax_rad+iMin_rad)/2;
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
        iCrop = int_rad(1:tidx_rad);
        tCrop_ax = t_ax(1:tidx_ax);
        fCrop_ax = FWHM_ax(1:tidx_ax);
        tCrop_rad = t_rad(1:tidx_rad);
    else
        dCrop = diams;
        iCrop = int_rad;
        tCrop_ax = t_ax;
        fCrop_ax = FWHM_ax;
        tCrop_rad = t_rad;
    end
    
    ud.rawDat(ii).cuttime = tCrop_rad;
    ud.rawDat(ii).cutdiams = dCrop;
    ud.rawDat(ii).cutint = iCrop;
    ud.rawDat(ii).cuttime_ax = tCrop_ax;
    ud.rawDat(ii).cutdiams_ax = fCrop_ax;

    if param.plot_raw
        plot(hh, tCrop_rad, dCrop, 'DisplayName', num2str(ii))
    end

end

ud.param = param;

if param.plot_raw
    set(gr, 'UserData', ud)
    
    if param.save_raw
        savefig(gr, fname)
    end
end
