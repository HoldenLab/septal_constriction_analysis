% Author: Kevin Whitley
% Date created: 190729

% This script asks the user for a directory containing a time-lapse image
% file (.tif) and a microbeJ track output (.mat). It then calls
% separate_microbej_tracks, which takes the centroid of each track over
% time and makes a new sub-image out of it. It then calls
% fit_septum_explicit to fit each one, and gets the diameters over time.

today = datestr(now, 'yymmdd');

path = uigetdir('D:\Data_Theia\');
dirIm = dir([path '\*.tif']);
dirTr = dir([path '\*microbej.mat']);

% fullstack = imreadstack([path '\' dirIm.name]);
% fullstack = im2double(fullstack);
% 
% tracks = load([path '\' dirTr.name]);

param.plot_im = 0; % plot septum images with fitted diameter
param.plot_raw = 0; % plot raw constriction vs. time data
param.plot_filt = 1; % filtered data (crap removed)
plot_t_con = 1; % histogram of constriction times

% params for separate_microbej_tracks
param.s_box = 34; % size of sub-images
param.n_frames_before = 5; % number of frames before first one to include

% params for fit_septum_explicit
param.psfFWHM = 250; % [nm]
param.pixSz = 65; % [nm/pix]
param.gm_model = 0; % use mixed Gaussian model (else - tilted circle model)
param.gridsp = 0.5;
param.fval_thresh = 0.0025;
if param.gm_model
    param.ccthresh = 0.965;
else
    param.ccthresh = 0.976;
end

orientation_cut = 0; % seems to not work so well. leave off for now.
param.interval = 2; % [min/frame]
param.fit_model = 1;

% params for fit_constriction_data
param.model = 'parabolic';
param.chi_thr = 1e4;

t_cpd = 40; % [min] time at which compound flushed in

% if param.plot_raw
%     figure('FileName', [path '/' today '_const_raw.fig'])
%     hh = gca;
%     hold on
%     box on
%     title(['Raw constriction traces for ' strrep(dirIm.name,'_','\_')])
%     xlabel(hh,'Time (min)')
%     ylabel(hh,'Septum width (nm)')
% end
% 
% if param.plot_filt
%     figure('FileName', [path '/' today '_const_filt.fig'])
%     gg = gca;
%     hold on
%     box on
%     legend('show')
%     title(['Filtered constriction traces for ' strrep(dirIm.name,'_','\_')])
%     xlabel(gg,'Time (min)')
%     ylabel(gg,'Septum width (nm)')
% end

[results, t_con] = fit_septum_explicit_batch_stripped(path, param);

% ntracks = length(tracks.Experiment.Lineage);

% subs=[]; norms=[]; fits=[]; results={}; rn=[]; se=[]; t_con=[];
% for ii = 1:ntracks
%     
%     %% Load image stack and pull out individual septa
%     
%     display(num2str(ii))
%     if ~isfield(tracks.Experiment.Lineage(ii).Trajectory(1).LIFESPAN,'f')
%         disp('trajectory removed in microbeJ.')
%         continue
%     elseif length(tracks.Experiment.Lineage(ii).Trajectory)>1
%         disp('microbeJ recorded division event. removing.')
%         continue
%     else
%         stack = separate_microbej_tracks(fullstack, tracks, ii, param.s_box, param.n_frames_before);
%     end
%     
%     %% Fit all frames to explicit septum model
%     
%     fits{ii} = fit_septum_explicit(stack, param.plot_im, param);
%     
%     start_frame = max([tracks.Experiment.Lineage(ii).FRAME.start - param.n_frames_before, 1]);
%     frames = (start_frame:tracks.Experiment.Lineage(ii).FRAME.end)';
%     frames = double(frames);
%     results{ii}.file = dirIm.name;
%     results{ii}.frames = frames;
%     results{ii}.fits = fits{ii};
%     
%     % Deprecated 190918
% %     if orientation_cut
% %         sept_or = atan((fits{ii}(:,4)-fits{ii}(:,2))./(fits{ii}(:,3)-fits{ii}(:,1))) * 180 / pi; % orientation of 'septa' determined by fits
% %         goodones = sept_or < -cell_or(1)+90+40 & sept_or > -cell_or(1)+90-40;
% %         frames(~goodones) = [];
% %         fits{ii}(~goodones,:) = []; % remove septa not orthogonal to cell axis
% %     end
%     
%     if isempty(frames)
%         continue
%     end
%     
%     if param.gm_model
%         subs = fits{ii}(:,1:2) - fits{ii}(:,3:4);
%         norms = sqrt(subs(:,1).^2 + subs(:,2).^2);
%     else
%         norms = fits{ii}(:,3)*param.interval;
%     end
%     results{ii}.diams = norms;
%     
%     frames(isnan(norms)) = [];
%     norms(isnan(norms)) = [];
%     
%     results{ii}.cutframes = frames;
%     results{ii}.cutdiams = norms;
%     
%     if plot_raw
%         plot(hh, frames*param.interval, norms*param.pixSz, 'DisplayName', num2str(ii))
%     end
%     
%     %% Fit trace of septum diameter vs. time to constriction model
%     
%     if fit_model
%         
%         ffit = frames(norms~=0) * param.interval; % [min]
%         lfit = norms(norms~=0) * param.pixSz; % [nm]
%         
%         % trim to remove ends of traces (typically not reliable)
%         indx = find(lfit<150, 1, 'first');
%         ffit((indx+1):end) = [];
%         lfit((indx+1):end) = [];
%         
%         if length(ffit)<6 || ~any(lfit<200) || ~any(lfit>900) % too short, or doesn't cover range of constriction - ignore
%             continue
%         end
% 
%         [results{ii}.cfit, rn(ii), se(ii,:), fcn] = fit_constriction_data(ffit, lfit, param.model);
%         
%         results{ii}.model = param.model;
%         results{ii}.fcn = fcn;
%         
%         results{ii}.chisq = rn(ii)/length(ffit);
%         
%         if ~isempty(ffit) && results{ii}.chisq<param.chi_thr
%             
%             t_con(ii) = ffit(end) - results{ii}.cfit(1);
%             
%             if plot_filt
%                 %         plot(gg, frames*2, norms*65)
%                 plot(gg, ffit, lfit, 'linew', 1, 'DisplayName', num2str(ii))
%                 tvf = ffit(1):0.1:ffit(end);
%                 plot(gg, tvf, fcn(results{ii}.cfit,tvf), 'k', 'DisplayName', '')
%                 
%                 results{ii}.time_filt = ffit;
%                 results{ii}.width_filt = lfit;
%             end
%         end
%     end
% end

%% Finish plotting results

% if param.plot_filt
% %     emptyCells = cellfun(@isempty,filtleg);
% %     filtleg = filtleg(~emptyCells);
% %     ggobj = ggobj(~emptyCells);
% %     legend(gg, ggobj, filtleg)
%     ud.analysis_date = today;
%     ud.im_name = dirIm.name;
%     ud.fname = dirTr.name;
%     ud.param = param;
%     ud.ntracks = ntracks;
%     
%     set(fhandle, 'UserData', ud)
% end

if plot_t_con
    figure('FileName', [path '/' today '_t_con.fig'])
    t_con(chisq>=param.chi_thr) = [];
    t_con(t_con<1) = [];
    hist(t_con)
    xlabel('t_{constriction} (min)')
    ylabel('Counts')
end
