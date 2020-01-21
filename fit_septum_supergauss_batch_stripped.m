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

function [ud, t_con] = fit_septum_supergauss_batch_stripped(fullstack, tracks, param)

if param.plot_raw
%     fname_base = extractBefore(param.im_file, '.ome');
%     fname_base = erase(fname_base, 'MMStack_');
%     gr = figure('FileName', [param.path '/' param.analysis_date '_' fname_base '_const_raw.fig']);
    gr = figure('FileName', [param.path '/' param.analysis_date]);
    hh = gca;
    hold on
    box on
    title(['Raw constriction traces for ' strrep(param.im_file,'_','\_')])
    xlabel(hh,'Time (min)')
    ylabel(hh,'Septum width (nm)')
end

if param.plot_filt
    fname_base = extractBefore(param.im_file, '.ome');
    fname_base = erase(fname_base, 'MMStack_');
    gf = figure('FileName', [param.path '/' param.analysis_date '_' fname_base '_const_filt.fig']);
    figure
    gg = gca;
    hold on
    box on
    title(['Filtered constriction traces for ' strrep(param.im_file,'_','\_')])
    xlabel(gg,'Time (min)')
    ylabel(gg,'Septum width (nm)')
end

subs=[]; diams=[]; fits=[]; rn=[]; t_con=[]; jj = 0; kk = 1;
for ii = 1:param.ntracks
    
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
    
    %% Fit all frames to explicit septum model
    
    [improfs, improfs_ax] = fit_septum_supergauss(stack, param.plot_im, param);
    
%     start_frame = tracks.Experiment.Lineage(ii).Trajectory.Bacteria(1).TRAJECTORY.frame.start + 1; % index in this file starts at zero, so add one
%     end_frame = tracks.Experiment.Lineage(ii).Trajectory.Bacteria(1).TRAJECTORY.frame.end + 1;
%     start_frame_new = max([start_frame - param.n_frames_before, 1]); % add a few extra frames before
%     frames = (start_frame_new:end_frame)';
    
    start_frame_new = max([imframes(1) - param.n_frames_before, 1]); % add a few extra frames before
    frames = ([(start_frame_new:start_frame_new+param.n_frames_before-1) imframes] + 1)';
    
    frames = double(frames);
    
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
    
    % need a way to delete all the crap after constriction finishes
    % Constriction end is just after intensity peak
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
    tidx = find(t==tConsEnd);
    tidx_ax = find(t_ax==tConsEnd_ax);
    
    tCrop = t(1:tidx);
    fCrop = FWHM(1:tidx);
    iCrop = intensity(1:tidx);
    tCrop_ax = t_ax(1:tidx_ax);
    fCrop_ax = FWHM_ax(1:tidx_ax);
    iCrop_ax = int_ax(1:tidx_ax);
    
    ud.imDat(ii).num = ii;
    ud.imDat(ii).line_profiles = improfs;
    ud.imDat(ii).frames = (imframes + 1)';
    ud.imDat(ii).centx = centx';
    ud.imDat(ii).centy = centy';
    
    % Deprecated 190918
    %     if orientation_cut
    %         sept_or = atan((fits{ii}(:,4)-fits{ii}(:,2))./(fits{ii}(:,3)-fits{ii}(:,1))) * 180 / pi; % orientation of 'septa' determined by fits
    %         goodones = sept_or < -cell_or(1)+90+40 & sept_or > -cell_or(1)+90-40;
    %         frames(~goodones) = [];
    %         fits{ii}(~goodones,:) = []; % remove septa not orthogonal to cell axis
    %     end
    
%     if isempty(time)
%         continue
%     end

    ud.rawDat(ii).num = ii;
    ud.rawDat(ii).width = FWHM;
    ud.rawDat(ii).time = t;
    ud.rawDat(ii).intensity = intensity; % added 200120
    ud.rawDat(ii).width_ax = FWHM_ax;
    ud.rawDat(ii).time_ax = t_ax;
    
%     time(isnan(diams)) = [];
%     diams(isnan(diams)) = [];
    
    ud.rawDat(ii).cuttime = tCrop;
    ud.rawDat(ii).cutdiams = fCrop;
    ud.rawDat(ii).cutint = iCrop;
    ud.rawDat(ii).cuttime_ax = tCrop_ax;
    ud.rawDat(ii).cutdiams_ax = fCrop_ax;
    
    if param.plot_raw
        plot(hh, tCrop, fCrop, 'DisplayName', num2str(ii))
    end
    
    %% Fit trace of septum diameter vs. time to constriction model
    
%     if param.fit_model
%         
%         tfit = time(diams~=0); % [min]
%         dfit = diams(diams~=0); % [nm]
%         
%         % trim to remove ends of traces (typically not reliable)
%         indx = find(dfit<200, 1, 'first');
%         tfit((indx+1):end) = [];
%         dfit((indx+1):end) = [];
%         
%         if length(tfit)<6 % too short - don't bother trying to fit
%             continue
%         end
% 
%         jj = jj + 1;
%         [cfit, rn(ii), se(ii,:), fcn] = fit_constriction_data(tfit, dfit, param.model);
%         
%         ud.fcn = fcn;
%         ud.fitDat(jj).num = ii;
%         ud.fitDat(jj).time = tfit;
%         ud.fitDat(jj).width = dfit;
%         ud.fitDat(jj).constFitVals = cfit;
%         
%         chisq = rn(ii)/length(tfit);
%         ud.fitDat(jj).constFitChiSq = rn(ii)/length(tfit);
%         
%         if ~isempty(tfit) && chisq<param.chi_thresh
%             
%             t_con(ii) = tfit(end) - cfit(1);
%             
%             if param.plot_filt
%                 plot(gg, tfit, dfit, 'linew', 1, 'DisplayName', num2str(ii))
%                 tvf = tfit(1):0.1:tfit(end);
%                 plot(gg, tvf, fcn(cfit,tvf), 'k', 'DisplayName', '')
%             end
%         end
%     end
end

ud.param = param;

set(gr, 'UserData', ud)

% if param.plot_filt
% 
%     set(gf, 'UserData', ud)
%     
%     tv = 0:1200;
%     plot(gg, ones(length(tv),1)*param.t_cpd, tv, 'r')
% end