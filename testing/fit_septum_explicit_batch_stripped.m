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

function [ud, t_con] = fit_septum_explicit_batch_stripped(fullstack, tracks, param)

if param.plot_raw
    figure('FileName', [param.path '/' param.analysis_date '_const_raw.fig'])
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
    gg = gca;
    hold on
    box on
    title(['Filtered constriction traces for ' strrep(param.im_file,'_','\_')])
    xlabel(gg,'Time (min)')
    ylabel(gg,'Septum width (nm)')
end

subs=[]; diams=[]; fits=[]; rn=[]; se=[]; t_con=[]; jj = 0; kk = 1;
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
    
    fits{ii} = fit_septum_explicit(stack, param.plot_im, param);
    
    start_frame = tracks.Experiment.Lineage(ii).Trajectory.Bacteria(1).TRAJECTORY.frame.start + 1; % index in this file starts at zero, so add one
    end_frame = tracks.Experiment.Lineage(ii).Trajectory.Bacteria(1).TRAJECTORY.frame.end + 1;
    start_frame_new = max([start_frame - param.n_frames_before, 1]); % add a few extra frames before
    frames = (start_frame_new:end_frame)';
    frames = double(frames);
    
    % add in excluded frames (that were out of focus or whatever)
    if ~isempty(param.exclude_frames) && any(frames(1)<param.exclude_frames) && any(frames(end)>param.exclude_frames)
        xx = param.exclude_frames;
        [exind, ~] = find(frames==xx);
        for kk = 1:length(exind)
            frames = [frames(1:exind(kk)-1); xx(kk); frames(exind(kk):end)+1];
            %         centx = [centx(1:exind-1) NaN centx(exind:end)];
            %         centy = [centy(1:exind-1) NaN centy(exind:end)];
            fits{ii} = [fits{ii}(1:exind(kk)-1,:); ones(1,4)*NaN; fits{ii}(exind(kk):end,:)];
        end
    end
    
    time = frames * param.interval; % [min]
    
    ud.imDat(ii).num = ii;
    ud.imDat(ii).imageFitVals = fits{ii};
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
    
    if isempty(time)
        continue
    end
    
    if param.gm_model
        subs = fits{ii}(:,1:2) - fits{ii}(:,3:4);
        diams_px = sqrt(subs(:,1).^2 + subs(:,2).^2); % [pixels]
    else
        diams_px = fits{ii}(:,3) * 2; % [pixels]
    end
    diams = diams_px * param.pixSz; % [nm]
    
    ud.rawDat(ii).num = ii;
    ud.rawDat(ii).width = diams;
    ud.rawDat(ii).time = time;
    
    time(isnan(diams)) = [];
    diams(isnan(diams)) = [];
    
    ud.rawDat(ii).cuttime = time;
    ud.rawDat(ii).cutdiams = diams;
    
    if param.plot_raw
        plot(hh, time, diams, 'DisplayName', num2str(ii))
    end
    
    %% Fit trace of septum diameter vs. time to constriction model
    
    if param.fit_model
        
        tfit = time(diams~=0); % [min]
        dfit = diams(diams~=0); % [nm]
        
        % trim to remove ends of traces (typically not reliable)
        indx = find(dfit<200, 1, 'first');
        tfit((indx+1):end) = [];
        dfit((indx+1):end) = [];
        
        if length(tfit)<6 % too short - don't bother trying to fit
            continue
        end

        jj = jj + 1;
        [cfit, rn(ii), se(ii,:), fcn] = fit_constriction_data(tfit, dfit, param.model);
        
        ud.fcn = fcn;
        ud.fitDat(jj).num = ii;
        ud.fitDat(jj).time = tfit;
        ud.fitDat(jj).width = dfit;
        ud.fitDat(jj).constFitVals = cfit;
        
        chisq = rn(ii)/length(tfit);
        ud.fitDat(jj).constFitChiSq = rn(ii)/length(tfit);
        
        if ~isempty(tfit) && chisq<param.chi_thresh
            
            t_con(ii) = tfit(end) - cfit(1);
            
            if param.plot_filt
                plot(gg, tfit, dfit, 'linew', 1, 'DisplayName', num2str(ii))
                tvf = tfit(1):0.1:tfit(end);
                plot(gg, tvf, fcn(cfit,tvf), 'k', 'DisplayName', '')
            end
        end
    end
end

ud.param = param;

if param.plot_filt

    set(gf, 'UserData', ud)
    
    tv = 0:1200;
    plot(gg, ones(length(tv),1)*param.t_cpd, tv, 'r')
end