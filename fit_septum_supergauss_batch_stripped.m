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

function [ud, t_con] = fit_septum_supergauss_batch_stripped(fullstack, tracks, param, tracknums)

if nargin < 4
    tracknums = 1:param.ntracks;
end

if param.plot_raw
    fname_base = extractBefore(param.im_file, '.ome');
    fname_base = erase(fname_base, 'MMStack_');
    gr = figure('FileName', [param.path '/' param.analysis_date '_' fname_base '_const_raw.fig']);
%     gr = figure('FileName', [param.path '/' param.analysis_date]);
    hh = gca;
    hold on
    box on
    title(['Raw constriction traces for ' strrep(param.im_file,'_','\_')])
    xlabel(hh,'Time (min)')
    ylabel(hh,'Septum width (nm)')
end

subs=[]; diams=[]; fits=[]; rn=[]; t_con=[]; jj = 0; kk = 1;
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
    
    [improfs, improfs_ax] = fit_septum_supergauss(stack, param.plot_im, param);
    
%     start_frame = tracks.Experiment.Lineage(ii).Trajectory.Bacteria(1).TRAJECTORY.frame.start + 1; % index in this file starts at zero, so add one
%     end_frame = tracks.Experiment.Lineage(ii).Trajectory.Bacteria(1).TRAJECTORY.frame.end + 1;
%     start_frame_new = max([start_frame - param.n_frames_before, 1]); % add a few extra frames before
%     frames = (start_frame_new:end_frame)';
    
    avail_frames_before = imframes(1); % how many frames before can we use?
    frames_before = min([param.n_frames_before avail_frames_before]);
    start_frame_new = max([imframes(1) - param.n_frames_before, 1]); % add a few extra frames before
    frames = ([(start_frame_new:start_frame_new+frames_before-1) imframes])'; % removed + 1 200124
%     frames = ([(start_frame_new:start_frame_new+frames_before-1) imframes] + 1)';
    
    frames = double(frames);
    
    ud.imDat(ii).num = ii;
    ud.imDat(ii).line_profiles = improfs;
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
%             fits{ii} = [fits{ii}(1:exind(kk)-1,:); ones(1,4)*NaN; fits{ii}(exind(kk):end,:)];
            improfs = [improfs(1:exind(kk)-1,:); ones(1,size(improfs,2))*NaN; improfs(exind(kk):end,:)];
            
            improfs_ax = [improfs_ax(1:exind(kk)-1,:); ones(1,size(improfs_ax,2))*NaN; improfs_ax(exind(kk):end,:)];
        end
    end
    
    time = frames * param.interval; % [min]
    
    %% Fit line profiles to generalized Gaussian function to get FWHM
    
    FWHM = []; intensity = []; widthResult = []; Rsq = []; se = []; FWHM_ax = []; se_ax = [];
    for jj = 1:size(improfs,1)
        linep = improfs(jj,~isnan(improfs(jj,:)));
        if isempty(linep) || length(linep)<5
            FWHM(jj) = NaN;
            intensity(jj) = NaN;
            widthResult(:,jj) = [NaN NaN NaN];
            se(:,jj) = ones(1,5)*NaN;
        else
            [FWHM(jj), intensity(jj), widthResult(:,jj), Rsq(jj), se(:,jj)] = fitProfile_public(linep,[],1);
        end
    end
    for kk = 1:size(improfs_ax,1)
        linep_ax = improfs_ax(kk,~isnan(improfs_ax(kk,:)));
        if isempty(linep_ax) || length(linep_ax)<5
            FWHM_ax(kk) = NaN;
            se_ax(:,kk) = ones(1,5)*NaN;
        else
            [FWHM_ax(kk), ~, ~, ~, se_ax(:,kk)] = fitProfile_public(linep_ax,[],1); % axial FWHM
        end
    end
    
    %% Cut results
    
    immin = time(1);
    immax = time(1)+size(improfs,1);

    MAXDIAM = 1500;
    MINDIAM = 50;
    FWHM = FWHM.*param.pixSz;
    FWHM_ax = FWHM_ax.*param.pixSz;
%     t = immin:immax;
    t = time;
    t_ax = time;
    int_ax = intensity;
    badIdx = FWHM>MAXDIAM | FWHM<MINDIAM | isnan(FWHM);
    badIdx_ax = FWHM_ax<MINDIAM | isnan(FWHM_ax);
    FWHM(badIdx) = [];
    intensity(badIdx) = [];
    t(badIdx) = [];
%     Rsq(badIdx) = [];
    se(:,badIdx) = [];
    FWHM_ax(badIdx_ax) = [];
    t_ax(badIdx_ax) = [];
    se_ax(:,badIdx_ax) = [];
    int_ax(badIdx_ax) = [];
    
%     badIdx2 = Rsq < 0.6;
%     FWHM(badIdx2) = [];
%     intensity(badIdx2) = [];
%     t(badIdx2) = [];
%     se(:, badIdx2) = [];
    
    badIdx3 = se(4,:) > 100;
    FWHM(badIdx3) = [];
    intensity(badIdx3) = [];
    t(badIdx3) = [];
    badIdx3_ax = se_ax(4,:) > 2;
    FWHM_ax(badIdx3_ax) = [];
    t_ax(badIdx3_ax) = [];
    int_ax(badIdx3_ax) = [];
    
    ud.rawDat(ii).num = ii;
    ud.rawDat(ii).width = FWHM;
    ud.rawDat(ii).time = t;
    ud.rawDat(ii).intensity = intensity; % added 200120
    ud.rawDat(ii).width_ax = FWHM_ax;
    ud.rawDat(ii).time_ax = t_ax;
    
    % need a way to delete all the crap after constriction finishes
    % Constriction end is just after intensity peak
    if param.intensity_cut
        [iMax, idxMax] = max(intensity);
        [iMax_ax, idxMax_ax] = max(int_ax);
        %     tMax= t(idxMax);
        
        [iMin, idx] = min(intensity(t>mean(t)));
        iMin_ax = min(int_ax(t_ax>mean(t_ax)));
        %     tMin=t(idx);
        iConsEnd = (iMax+iMin)/2;
        tCrop = t(idxMax:end);
        iCrop = intensity(idxMax:end);
        iConsEnd_ax = (iMax_ax+iMin_ax)/2;
        tCrop_ax = t_ax(idxMax_ax:end);
        iCrop_ax = int_ax(idxMax_ax:end);
        
        iDiff = iCrop-iConsEnd;
        iConsEnd = iCrop(find(iDiff<0,1)-1);
        tConsEnd = tCrop(find(iDiff<0,1)-1);
        iDiff_ax = iCrop_ax - iConsEnd_ax;
        tConsEnd_ax = tCrop_ax(find(iDiff_ax<0,1)-1);
        
        if isempty(tConsEnd)
            continue
        end
        if isempty(tConsEnd_ax)
            continue
        end
        tidx = find(t==tConsEnd);
        tidx_ax = find(t_ax==tConsEnd_ax);
        
        tCrop = t(1:tidx);
        fCrop = FWHM(1:tidx);
        iCrop = intensity(1:tidx);
        tCrop_ax = t_ax(1:tidx_ax);
        fCrop_ax = FWHM_ax(1:tidx_ax);
        iCrop_ax = int_ax(1:tidx_ax);
    else
        tCrop = t;
        fCrop = FWHM;
        iCrop = intensity;
        tCrop_ax = t_ax;
        fCrop_ax = FWHM_ax;
    end
    
    ud.rawDat(ii).cuttime = tCrop;
    ud.rawDat(ii).cutdiams = fCrop;
    ud.rawDat(ii).cutint = iCrop;
    ud.rawDat(ii).cuttime_ax = tCrop_ax;
    ud.rawDat(ii).cutdiams_ax = fCrop_ax;

%     time(isnan(diams)) = [];
%     diams(isnan(diams)) = [];
    
    if param.plot_raw
        plot(hh, tCrop, fCrop, 'DisplayName', num2str(ii))
    end

end

ud.param = param;

set(gr, 'UserData', ud)
