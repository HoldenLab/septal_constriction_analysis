% Author: Kevin Whitley
% Date created: 191106

% This script is used for fitting and plotting various data from
% constriction experiments. Input is a structure variable called alldat2,
% which must contain the field rawDat.

filter_dat = 1; % filter data based on fitting error, etc.
fit_dat = 1; % process and fit data to selected model for constriction
model = 'parabolic'; % constriction model to use
Rsq_thresh = 0.8; % R^2 threshold for data to be plotted. also threshold for fitted values to be included in teff violin plots.
d0 = 1100; % fix d0 for calculating teff. updated 200625.

plot_raw_dat = 0; % plot raw data traces
plot_filtered_dat = 0; % plot filtered data traces
plot_fitted_dat = 1; % plot processed traces with fitted model
plot_teff = 0; % make violin plots of effective constriction time (teff)
plot_d_vs_alph = 0; % plot diameter at time of drug treatment vs. constriction rate
plot_perturbed = 0; % plot only traces where data was recorded DURING drug treatment
plot_intensity = 0; % plot septal intensities over time
plot_width_axial = 0; % plot axial widths (thickness) over time
plot_d0 = 0; % plot fitted d0 (pretty useless)

%% Load variable dat, which has all raw data from alldat2
dat=[]; undat=[];
for ii = 1:length(alldat2)
    for jj = 1:length(alldat2(ii).rawDat)
        
        rawdat = alldat2(ii).rawDat(jj);
%         unfiltdat = alldat2(ii).datUnfilt(jj);
        
        unfiltdat.num = alldat2(ii).rawDat(jj).num;
        
%         rawdat.param = alldat2(ii).param;
%         rawdat.im_date = alldat2(ii).param.im_file(1:6);
        unfiltdat.param = alldat2(ii).param;
        unfiltdat.im_date = alldat2(ii).param.im_file(1:6);
        
        posind = strfind(alldat2(ii).param.im_file, 'Pos');
%         rawdat.pos = alldat2(ii).param.im_file(posind:posind+3);
        unfiltdat.pos = alldat2(ii).param.im_file(posind:posind+3);
        
        dat = [dat; rawdat];
%         undat = [undat; unfiltdat];
%         dat = [dat; unfiltdat];
    end
end

%% Plot all raw data traces (no filters, no fits)
if plot_raw_dat
    figure('FileName', [path '/' today '_raw_dat.fig'])
    hold on
    box on
    
    for ii = 1:length(dat)
        plot(dat(ii).cuttime, dat(ii).cutdiams, ...
            'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num)])
    end
    tv = 0:1200;
    plot(ones(length(tv),1)*dat(ii).param.t_cpd, tv, 'r')
    
    xlabel('Time (min)')
    ylabel('Diameter (nm)')
    title('Raw data')
end

%% Filter data
if filter_dat
    
    % thresholds
%     if isfield(undat(1),'linecents') % for fitted line centers
%         allcents = cat(1,undat(:).linecents);
    if isfield(dat(1),'linecents') % for fitted line centers
        allcents = cat(1,dat(:).linecents);
        centmin = nanmean(allcents)-2*nanstd(allcents);
        centmax = nanmean(allcents)+2*nanstd(allcents);
        %         centmin = 15;
        %         centmax = 23.5;
        %         centmin = 11;
        %         centmax = 13.1;
    end
    
    minthresh_ax = 0; % min thickness (unused)
    maxthresh_ax = 1300; % max thickness (unused)
    %         errthresh_ax = 0.7; % for standard error in one fitted parameter
    errthresh_ax = 0.90; % for R^2
            errthresh_rad = 0.9; % for residual norm
    %         errthresh_rad = 1.5; % for residual norm
%     errthresh_rad = 0.90; % for R^2
    
    for ii = 1:length(dat)
%         tdat = undat(ii).time(:);
%         axdat = undat(ii).FWHM_ax(:);
%         edat_ax = undat(ii).se_ax(:);
%         raddat = undat(ii).diam(:);
%         edat_rad = undat(ii).fiterrs(:);
        tdat = dat(ii).time(:);
        axdat = dat(ii).FWHM_ax(:);
        edat_ax = dat(ii).se_ax(:);
        raddat = dat(ii).diam(:);
        edat_rad = dat(ii).fiterrs(:);
        
        % cut on fitted center of line profile (if center close to edge of
        % image it often gives huge diameters with artificially small
        % errors).
%         if isfield(undat(ii),'linecents')
%             centdat = undat(ii).linecents(:);
        if isfield(dat(ii),'linecents')
            centdat = dat(ii).linecents(:);
            
            centind = centdat>centmin & centdat<centmax;
            tdat = tdat(centind);
            edat_ax = edat_ax(centind);
            axdat = axdat(centind);
            raddat = raddat(centind);
            edat_rad = edat_rad(centind);
        end

        % cut axial and radial times separately
        tdat_rad = tdat;
        tdat_ax = tdat;

        % cut thickness on max/min values - disabled, usually not needed
        if 0
            tdat_ax = tdat_ax(axdat>minthresh_ax & axdat<maxthresh_ax);
            edat_ax = edat_ax(axdat>minthresh_ax & axdat<maxthresh_ax);
            axdat = axdat(axdat>minthresh_ax & axdat<maxthresh_ax);
        end
        
        % cut on error in thickness fit. data in different formats -
        % older format was standard error in one parameter fit, new format
        % is R^2
%         tdat_ax = tdat_ax(edat_ax<errthresh_ax); % se
%         axdat = axdat(edat_ax<errthresh_ax); % se
        tdat_ax = tdat_ax(edat_ax>errthresh_ax); % R^2
        axdat = axdat(edat_ax>errthresh_ax); % R^2
        
        % cut on error in diameter fit. data in different formats - older
        % format was residual norm, newer is R^2
        tdat_rad = tdat_rad(edat_rad<errthresh_rad); % resnorm
        raddat = raddat(edat_rad<errthresh_rad); % resnorm
%         tdat_rad = tdat_rad(edat_rad>errthresh_rad); % R^2
%         raddat = raddat(edat_rad>errthresh_rad); % R^2
 
        % cut flat ends of traces (lingering intensity, not realy
        % constriction) using first derivative
        d_s = smooth(raddat,5); % smooth factor was 10. now 5 (200515).
        t_s = smooth(tdat_rad,5);
        dDdt = (d_s(2:end)-d_s(1:end-1)) ./ (t_s(2:end)-t_s(1:end-1)); % first derivative of diameters
        dDdt_s = smooth(dDdt,3); % smoothed first derivative. factor was 10. now 3 (200515).
        dDdt_thresh = 0; % threshold for first derivative
%         halfind = floor(length(dDdt_s)/2) + 1; % halfway mark of trace (can't remove any flat parts at beginning). +1 needed to handle super short traces.
        quartind = floor(length(dDdt_s)/4) + 1; % quarterway mark of trace (can't remove any flat parts at beginning). +1 needed to handle super short traces.
%         indx = min([find(dDdt_s(halfind:end)>=dDdt_thresh, 1, 'first')+halfind-1 length(tfit)]);
        indx = min([find(dDdt_s(quartind:end)>=dDdt_thresh, 1, 'first')+quartind-1 length(tdat_rad)]);
        
        % second derivative (deprecated)
        %         bins2 = bins(1:end-1) + (bins(2:end)-bins(1:end-1))/2;
        %         w2d = (w1d_s(2:end)-w1d_s(1:end-1))./bins2; % second derivative
        %         t1 = t_s(1:end-1) + (t_s(2:end)-t_s(1:end-1))/2;
        %         t2 = t1(1:end-1) + (t1(2:end)-t1(1:end-1))/2;
        %         w2d_thresh = 10;
        %         w2d_c = w2d(2:end-1);
        %         t2_c = t2(2:end-1);
        %         indx = min([find(w2d_c>=w2d_thresh, 1, 'first') length(ffit)]);
        
        tdat_rad_c = tdat_rad(1:indx); % cut times (flat ends removed)
        raddat_c = raddat(1:indx); % cut diameters
        if ~isempty(tdat_rad_c)
            tdat_ax_c = tdat_ax(tdat_ax<tdat_rad_c(end));
            axdat_c = axdat(tdat_ax<tdat_rad_c(end));
        end
          
        dat(ii).cuttime_ax = tdat_ax_c;
        dat(ii).cutdiams_ax = axdat_c;
        dat(ii).cuttime = tdat_rad_c;
        dat(ii).cutdiams = raddat_c;
    end
else
    for ii = 1:length(dat)
        dat(ii).cuttime_ax = dat(ii).time(~isnan(dat(ii).FWHM_ax));
        dat(ii).cutdiams_ax = dat(ii).FWHM_ax(~isnan(dat(ii).FWHM_ax));
        dat(ii).cuttime = dat(ii).time(~isnan(dat(ii).diam));
        dat(ii).cutdiams = dat(ii).diam(~isnan(dat(ii).diam));
    end
end

%% Plot filtered data (no fits)
if plot_filtered_dat
    figure('FileName', [path '/' today '_filtered_dat.fig'])
    
    for ii = 1:length(dat)

        ax1 = subplot(211);
        hold(ax1, 'on')
        plot(ax1, dat(ii).cuttime, dat(ii).cutdiams, 'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num)])
        ylim([200 1200])
        
        ax2 = subplot(212);
        hold(ax2,'on')
        plot(dat(ii).cuttime_ax, dat(ii).cutdiams_ax, 'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num)])
        
        % partition thickness traces, plot separately
        %         [cond, switchPt] = detectZRcondensation(dat(ii).cutdiams_ax, 'MinSpatialDist', 50, 'Statistic', 'linear');
        %         hold(ax2,'off')
        %         plot(ax2,dat(ii).cuttime_ax(cond==1), dat(ii).cutdiams_ax(cond==1), 'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num)], 'Color', 'k')
        %         plot(ax2,dat(ii).cuttime_ax(cond==2), dat(ii).cutdiams_ax(cond==2), 'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num)], 'Color', 'r')
        %         hold(ax2,'off')
        ylim([200 1000])
    end
    linkaxes([ax1 ax2],'x')
    tv = 0:1000;
    plot(ones(length(tv),1)*dat(ii).param.t_cpd, tv, 'r')
    
    title(ax1, 'Filtered data')
    ylabel(ax1, 'Diameter (nm)')
    xlabel(ax2, 'Time (min)')
    ylabel(ax2, 'Thickness (nm)')
end

%% Process and fit data to selected model for constriction
if fit_dat
    
    if plot_fitted_dat
        figure('FileName', [path '/' today '_fit_dat.fig'])
        ax1 = subplot(211);
        hold(ax1, 'on')
        box(ax1, 'on')
        ax2 = subplot(212);
        hold(ax2, 'on')
        box(ax2, 'on')
    end
    
    diffp=[]; diffl=[]; t_switch=[]; dax_cond=[]; dax_decond=[]; dax_precons=[]; dax_postcons=[];
    for ii = 1:length(dat)
        
        % 1. prepare data for fitting
        
        % make sure these are both column vectors
        tfit = dat(ii).cuttime(:);
        dfit = dat(ii).cutdiams(:);
        tax_c = dat(ii).cuttime_ax(:);
        dax_c = dat(ii).cutdiams_ax(:);

        % moved section to filter area 201020
        if 0
            % trim to remove flat ends of traces based on first derivative
            % (this is not constriction).
            d_s = smooth(dfit,5); % smooth factor was 10. now 5 (200515).
            t_s = smooth(tfit,5);
            dDdt = (d_s(2:end)-d_s(1:end-1)) ./ (t_s(2:end)-t_s(1:end-1)); % first derivative of diameters
            dDdt_s = smooth(dDdt,3); % smoothed first derivative. factor was 10. now 3 (200515).
            dDdt_thresh = 0; % threshold for first derivative
            %         halfind = floor(length(dDdt_s)/2) + 1; % halfway mark of trace (can't remove any flat parts at beginning). +1 needed to handle super short traces.
            quartind = floor(length(dDdt_s)/4) + 1; % quarterway mark of trace (can't remove any flat parts at beginning). +1 needed to handle super short traces.
            %         indx = min([find(dDdt_s(halfind:end)>=dDdt_thresh, 1, 'first')+halfind-1 length(tfit)]);
            indx = min([find(dDdt_s(quartind:end)>=dDdt_thresh, 1, 'first')+quartind-1 length(tfit)]);
            
            % second derivative (deprecated)
            %         bins2 = bins(1:end-1) + (bins(2:end)-bins(1:end-1))/2;
            %         w2d = (w1d_s(2:end)-w1d_s(1:end-1))./bins2; % second derivative
            %         t1 = t_s(1:end-1) + (t_s(2:end)-t_s(1:end-1))/2;
            %         t2 = t1(1:end-1) + (t1(2:end)-t1(1:end-1))/2;
            %         w2d_thresh = 10;
            %         w2d_c = w2d(2:end-1);
            %         t2_c = t2(2:end-1);
            %         indx = min([find(w2d_c>=w2d_thresh, 1, 'first') length(ffit)]);
            
            tfit_c = tfit(1:indx); % cut times (flat ends removed)
            dfit_c = dfit(1:indx); % cut diameters
            if ~isempty(tfit_c)
                tax_c = dat(ii).cuttime_ax(dat(ii).cuttime_ax<tfit_c(end));
                dax_c = dat(ii).cutdiams_ax(dat(ii).cuttime_ax<tfit_c(end));
            end
        end
        
        % partition thickness traces into decondensed and condensed
        % using state change
        [cond, switchPt, muse] = detectZRcondensation(dax_c, 'MinSpatialDist', 50, 'Statistic', 'std'); % from point change
        dax_cond = [dax_cond; dax_c(cond==1)];
        dax_decond = [dax_decond; dax_c(cond==2)];
        
        % using threshold
        thick_cut = 422; % [nm] mean + 3*sigma for 'condensed' state of wt FtsZ-GFP
        dax_c_sm = smooth(dax_c, 3);
        inds = dax_c_sm < thick_cut;
        trans = inds(1:end-1) - inds(2:end);
        transind = [0; find(trans==-1); find(trans==1); length(inds)];
        transind = sort(transind);

        % Split data into pre- and post-treatment
        % exception: if param.t_cpd=0, assume all data is untreated (or
        % 'pre-treated')
        if dat(ii).param.t_cpd == 0
            % modified 201020
%             dat(ii).preTime = tfit;
%             dat(ii).preDiam = dfit;
%             dat(ii).postTime = [];
%             dat(ii).postDiam = [];
            
%             dat(ii).preTimeCut = tfit_c;
%             dat(ii).preDiamCut = dfit_c;
%             dat(ii).postTimeCut = [];
%             dat(ii).postDiamCut = [];
            
            dat(ii).preTimeCut = tfit;
            dat(ii).preDiamCut = dfit;
            dat(ii).postTimeCut = [];
            dat(ii).postDiamCut = [];
        else
            % modified 201020
%             dat(ii).preTime = tfit(tfit<dat(ii).param.t_cpd);
%             dat(ii).preDiam = dfit(tfit<dat(ii).param.t_cpd);
%             dat(ii).postTime = tfit(tfit>=dat(ii).param.t_cpd);
%             dat(ii).postDiam = dfit(tfit>=dat(ii).param.t_cpd);
%             
%             dat(ii).preTimeCut = tfit_c(tfit_c<dat(ii).param.t_cpd);
%             dat(ii).preDiamCut = dfit_c(tfit_c<dat(ii).param.t_cpd);
%             dat(ii).postTimeCut = tfit_c(tfit_c>=dat(ii).param.t_cpd);
%             dat(ii).postDiamCut = dfit_c(tfit_c>=dat(ii).param.t_cpd);
            
            dat(ii).preTimeCut = tfit(tfit<dat(ii).param.t_cpd);
            dat(ii).preDiamCut = dfit(tfit<dat(ii).param.t_cpd);
            dat(ii).postTimeCut = tfit(tfit>=dat(ii).param.t_cpd);
            dat(ii).postDiamCut = dfit(tfit>=dat(ii).param.t_cpd);
        end
        
        % 2. fit the data
        
        % fit pre-treatment data
        if length(dat(ii).preTimeCut) >= 5 % needs to have enough points for decent fit
            [cfit_pre, rn_pre, res, fcnp] = fit_constriction_data(dat(ii).preTimeCut, dat(ii).preDiamCut, model);
            [cfit_prel, rn_prel, resl, fcnl] = fit_constriction_data(dat(ii).preTimeCut, dat(ii).preDiamCut, 'linear');
            
            % calculate goodness-of-fit metrics
            sst_pre = sum((dat(ii).preDiamCut - mean(dat(ii).preDiamCut)).^2); % sum of squares total
            Rsq_pre = 1 - rn_pre/sst_pre; % R^2
            Xsq_pre = rn_pre ./ length(dat(ii).preDiamCut); % chi^2 (kind of)
            
            % calculate goodness-of-fit for data only after constriction
            % begins (flat part at beginning really shouldn't contribute to
            % fitting error).
            onlyconpre = dat(ii).preDiamCut(dat(ii).preTimeCut >= cfit_pre(1));
            t_onlyconpre = dat(ii).preTimeCut(dat(ii).preTimeCut >= cfit_pre(1));
            n_pre = length(onlyconpre);
            rn_onlyconpre = sum(res(dat(ii).preTimeCut >= cfit_pre(1)).^2);
            sst_pre_con = sum((onlyconpre - mean(onlyconpre)).^2);
            Rsq_pre_con = 1 - rn_onlyconpre/sst_pre_con;
            Xsq_pre_con = sum((onlyconpre-fcnp(cfit_pre,t_onlyconpre)).^2) / length(onlyconpre);
        else
            cfit_pre = [];
            fcnp = [];
            Rsq_pre = NaN;
            Rsq_pre_con = NaN;
            n_pre = NaN;
            Xsq_pre = NaN;
            Xsq_pre_con = NaN;
        end
        
        dat(ii).fitDat.preFitVals = cfit_pre;
        dat(ii).fitDat.preRsq = Rsq_pre;
        dat(ii).fitDat.preRsq_con = Rsq_pre_con;
        dat(ii).fitDat.n_pre = n_pre;
        dat(ii).fitDat.XSq = Xsq_pre;
        dat(ii).fitDat.XSq_pre_con = Xsq_pre_con;
        
        % fit post-treatment data
        if length(dat(ii).postTimeCut) >= 5 && any(dat(ii).postTimeCut <= dat(ii).param.t_cpd+1) % bit shoddy, but should work for now
            [cfit_post, rn_post, res, fcnp] = fit_constriction_data(dat(ii).postTimeCut, dat(ii).postDiamCut, model);
            
            % calculate goodness-of-fit metrics
            sst_post = sum((dat(ii).postDiamCut - mean(dat(ii).postDiamCut)).^2); % sum of squares total
            Rsq_post = 1 - rn_post/sst_post; % R^2
            
            % calculate goodness-of-fit for data only after constriction
            % begins (flat part at beginning really shouldn't contribute to
            % fitting error).
            onlycon = dat(ii).postDiamCut(dat(ii).postTimeCut >= cfit_post(1));
            n_post = length(onlycon);
            rn_onlycon = sum(res(dat(ii).postTimeCut >= cfit_post(1)).^2);
            sst_post_con = sum((onlycon - mean(onlycon)).^2);
            Rsq_post_con = 1 - rn_onlycon/sst_post_con;
        else
            cfit_post = [];
            Rsq_post = NaN;
            Rsq_post_con = NaN;
            n_post = NaN;
        end

        dat(ii).fcn = fcnp;
        dat(ii).fitDat.num = ii;
        
        dat(ii).fitDat.postFitVals = cfit_post;
        dat(ii).fitDat.postRsq = Rsq_post;
        dat(ii).fitDat.postRsq_con = Rsq_post_con;
        dat(ii).fitDat.n_post = n_post;

        % 3. plot the data with fits
        
        if plot_fitted_dat
            
            % plot data pre-treatment
            if Rsq_pre_con > Rsq_thresh && n_pre>5
                
                shift = dat(ii).param.t_cpd; % define t=0 as treatment time
%                 shift = cfit_pre(1); % define t=0 as constriction start time
                
                % hold(ax1,'off')
                % hold(ax2,'off')
                
                % plot diameters, shifted by either treatment time or
                % constriction start time
                plot(ax1, dat(ii).preTimeCut-shift, dat(ii).preDiamCut, 'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num) '\_pre'], 'linew', 1)
                
                % plot constriction model on top of data
                % hold(ax1,'on')
%                 tvf = dat(ii).preTime(1):0.1:dat(ii).preTime(end);
%                 plot(ax1, tvf-shift, fcnp(cfit_pre,tvf), 'k', 'DisplayName', '')
%                 plot(ax1, tvf-shift, fcnl(cfit_prel,tvf), 'k', 'DisplayName', '')

                % segment and plot thickness traces
                seg_t={}; seg_d={}; dur=[];
                for kk = 1:length(transind)-1
                    seg_t{kk} = tax_c(transind(kk)+1:transind(kk+1));
                    seg_d{kk} = dax_c(transind(kk)+1:transind(kk+1));
                    if ~isempty(tax_c) && transind(kk)+1<tax_c(end) && inds(transind(kk)+1)==0 % decondensed
                        clr = 'r';
                    elseif ~isempty(tax_c) % condensed
                        clr = 'k';
                        
                        % find every consecutive stretch of 'condensed'
                        % state prior to constriction, measure length of
                        % segment
                        precons = seg_t{kk}-shift;
                        precons = precons(precons<=0);
                        dur(kk) = length(precons);
%                         dur(kk) = length(seg_t{kk});
                    end
%                     plot(ax2, seg_t{kk}-shift, seg_d{kk}, 'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num) '\_pre'], 'linew', 1, 'Color', clr)
%                     hold(ax2,'on')
                end
                
                dat(ii).cond_dur = max([0 dur]);

                allprecons = dat(ii).preTimeCut-shift;
                allprecons = allprecons(allprecons<=0); % need to be sure there are actually 5 points there to compare
                if length(allprecons)>=1 && ~isempty(dat(ii).cond_dur) && dat(ii).cond_dur>=1 % condensation for 5 points pre-constriction
                    dat(ii).didcond = 1;
                elseif length(allprecons)>=1 && ~isempty(dat(ii).cond_dur) && dat(ii).cond_dur<1 % no condensation for 5 points pre-constriction (but there were 5 data points)
                    dat(ii).didcond = 0;
                end
                
                dax_cond = dax_c(dax_c_sm<thick_cut);
                dax_decond = dax_c(dax_c_sm>=thick_cut);
                
                dax_precons = [dax_precons; dax_c(tax_c<cfit_pre(1))];
                dax_postcons = [dax_postcons; dax_c(tax_c>=cfit_pre(1))];
                
%                 plot(ax2, tax_c(cond==1)-shift, dax_c(cond==1), 'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num) '\_pre'], 'linew', 1, 'Color', 'k')
%                 plot(ax2, tax_cond-shift, dax_cond, 'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num) '\_pre'], 'linew', 1, 'Color', 'k')
%                 hold(ax2,'on')
%                 plot(ax2, tax_c(cond==2)-shift, dax_c(cond==2), 'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num) '\_pre'], 'linew', 1, 'Color', 'r')
%                 plot(ax2, tax_decond-shift, dax_decond, 'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num) '\_pre'], 'linew', 1, 'Color', 'r')
                plot(ax2, tax_c-shift, dax_c, 'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num) '\_pre'], 'linew', 1)
                if switchPt>0
                    t_switch(ii) = tax_c(switchPt)-cfit_pre(1);
                end

                % calculate difference between 2B disappearance time and
                % predicted end of constriction. Since d0 is not constant
                % for each trace, can't use simple formula for end time of
                % constriction.
%                 gone2b = dat(ii).cuttime(end);
%                 tvf2 = dat(ii).preTime(1):0.01:dat(ii).preTime(end)+50;
%                 pararray = fcnp(cfit_pre, tvf2);
%                 teffp = d0^2 ./ cfit_pre(2);
%                 parzero = find(~pararray, 1, 'first');
%                 teffp = tvf2(parzero); % can't use fzero because lot of zeros
%                 teffl = -d0 / cfit_prel(2);
%                 teffl = fzero(@(x)fcnl(cfit_prel,x), cfit_pre(1));
                
%                 diffp = [diffp; teffp-gone2b];
%                 diffl = [diffl; teffl-gone2b];
            end
            
            % plot data post-treatment (if exists)
            if Rsq_post_con > Rsq_thresh && n_post>5
                
                plot(ax1, dat(ii).postTimeCut-shift, dat(ii).postDiamCut, 'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num) '\_post'], 'linew', 1)

                % plot constriction model on top of data
%                 tvf = dat(ii).postTimeCut(1):0.1:dat(ii).postTimeCut(end);
%                 plot(ax1, tvf-shift, fcnp(cfit_post,tvf), 'r', 'DisplayName', '')
%                 teff = d0^2 ./ cfit_post(2);
                
            end
            
        end
    end
    linkaxes([ax1 ax2], 'x')
end

%% Make violin plots (with DABEST if possible) of effective constriction times
if plot_teff

    teff_pre = []; teff_post = []; cat_pre = {}; cat_post = {};
    for ii = 1:length(dat)

        if ~isempty(dat(ii).fitDat) && dat(ii).fitDat.preRsq_con > Rsq_thresh && dat(ii).fitDat.n_pre>5
%             teff_pre(ii) = dat(ii).fitDat.preFitVals(3)^2 ./ dat(ii).fitDat.preFitVals(2);
            teff_pre = [teff_pre d0^2 ./ dat(ii).fitDat.preFitVals(2)];
            cat_pre = [cat_pre 'Untreated'];
%         else
%             teff_pre(ii) = [];
        end
        
        if ~isempty(dat(ii).fitDat) && dat(ii).fitDat.postRsq_con > Rsq_thresh && dat(ii).fitDat.n_post>5
%             teff_post(ii) = dat(ii).fitDat.postFitVals(3)^2 ./ dat(ii).fitDat.postFitVals(2);
            teff_post = [teff_post d0^2 ./ dat(ii).fitDat.postFitVals(2)];
            cat_post = [cat_post 'Treated'];
%         else
%             teff_post(ii) = [];
        end
    end
    
    teffs = [teff_pre teff_post]';
    cats = [cat_pre cat_post]';
    
%     teff_pre = [teff_pre; ones(length(teff_post)-length(teff_pre),1)*NaN];
%     teff_post = [teff_post; ones(length(teff_pre)-length(teff_post),1)*NaN];
    
    if any(~isnan(teff_post))
        violinplusdabest(teffs, cats, {'Untreated', 'Treated'})
    else
        figure('FileName', [path '/' today '_teff.fig'])
        
        violinplot(teff_pre)
        
        ylabel('t_{eff} (min)')
    end
end

%% Plot diameter of ring at time of treatment vs. constriction rate
if plot_d_vs_alph
    figure('FileName', [path '/' today '_alpha.fig'])
    hold on
    box on
    
    d_cpd = []; alph_cpd = [];
    for ii = 1:length(dat)
        if ~isempty(dat(ii).fitDat) && any(dat(ii).fitDat.time<dat(ii).param.t_cpd) && any(dat(ii).fitDat.time>dat(ii).param.t_cpd) && dat(ii).fitDat.postRsq > Rsq_thresh
            rng = 1; % [min]
            ind = find(dat(ii).fitDat.time>=dat(ii).param.t_cpd-rng & dat(ii).fitDat.time<=dat(ii).param.t_cpd+rng);
            d_all = dat(ii).fitDat.width(ind);
            d_cpd = [d_cpd; mean(d_all)];
            teff = dat(ii).fitDat.constFitVals(3)^2 ./ dat(ii).fitDat.constFitVals(2);
            alph_cpd = [alph_cpd; dat(ii).fitDat.constFitVals(2) teff];
        end
    end
    
    d_cpd(d_cpd==0) = [];
    alph_cpd(alph_cpd==0) = [];
    
    plot(d_cpd, alph_cpd(:,1), 'ok')
    xlabel('d_{constriction} (nm)')
    ylabel('\alpha (nm^2/min)')
end

%% Plot only traces where data was recorded DURING drug treatment
if plot_perturbed
    figure('FileName', [path '/' today '_perturbed.fig'])
    hold on
    box on
    
    for ii = 1:length(dat)
        %         if ~isempty(dat(ii).fitDat) && any(dat(ii).fitDat.time<dat(ii).param.t_cpd) && any(dat(ii).fitDat.time>dat(ii).param.t_cpd)
        %         if ~isempty(dat(ii).fitDat) && dat(ii).fitDat.time(end)<dat(ii).param.t_cpd
        if any(dat(ii).cuttime<dat(ii).param.t_cpd) && any(dat(ii).cuttime>dat(ii).param.t_cpd)
            %             if dat(ii).fitDat.constFitChiSq < chiSq_thresh
            plot(dat(ii).cuttime, dat(ii).cutdiams, 'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num)])
            %                 plot(dat(ii).fitDat.time-dat(ii).param.t_cpd, dat(ii).fitDat.width, 'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num)], 'linew', 1)
            
            %                 tvf = dat(ii).fitDat.time(1):0.1:dat(ii).fitDat.time(end);
            %                 plot(tvf-dat(ii).param.t_cpd, dat(ii).fcn(dat(ii).fitDat.constFitVals,tvf), 'k', 'DisplayName', '')
            %             end
        end
    end
end

%% Plot septal intensites vs. time
if plot_intensity
    figure('FileName', [path '/' today '_pre.fig'])
    subplot(212)
    h1 = gca;
    hold on
    subplot(211)
    g1 = gca;
    hold on
    
    figure('FileName', [path '/' today '_post.fig'])
    subplot(212)
    h2 = gca;
    hold on
    subplot(211)
    g2 = gca;
    hold on
    
    for ii = 1:length(dat)
        
        if ~any(dat(ii).cuttime>dat(ii).param.t_cpd) && ~isempty(dat(ii).fitDat.preFitVals) && dat(ii).fitDat.preRsq_con > Rsq_thresh && dat(ii).fitDat.n_pre>5 % only before
%         if any(dat(ii).cuttime<dat(ii).param.t_cpd) && any(dat(ii).cuttime>dat(ii).param.t_cpd) % only perturbed
%         if isfield(dat(ii).fitDat,'preFitVals') && ~isempty(dat(ii).fitDat.preFitVals)

            plot(g1, dat(ii).cuttime-dat(ii).fitDat.preFitVals(1), dat(ii).cutdiams, ...
                'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num)]) % post aligned by treatment time

            plot(h1, dat(ii).cuttime-dat(ii).fitDat.preFitVals(1), dat(ii).cutint, ...
                'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num)]) % pre aligned by constriction time
        end
        
%         if any(dat(ii).cuttime>dat(ii).param.t_cpd) && ~isempty(dat(ii).fitDat.postFitVals) && dat(ii).fitDat.postRsq_con > Rsq_thresh && dat(ii).fitDat.n_post>5 % only after
%         if any(dat(ii).cuttime>dat(ii).param.t_cpd) && any(dat(ii).cuttime>dat(ii).param.t_cpd) && isempty(dat(ii).fitDat.postFitVals) && ~any(dat(ii).cutdiams<700)% only after, raw data
        if any(dat(ii).cuttime==dat(ii).param.t_cpd) && dat(ii).cutdiams(dat(ii).cuttime==dat(ii).param.t_cpd)>900 % only after, mature rings only
            
%             plot(dat(ii).cuttime-dat(ii).fitDat.postFitVals(1), dat(ii).cutint, ...
%                 'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num)]) % post aligned by constriction time

            plot(g2, dat(ii).cuttime-dat(ii).param.t_cpd, dat(ii).cutdiams, ...
                'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num)]) % post aligned by treatment time
            
            plot(h2, dat(ii).cuttime-dat(ii).param.t_cpd, dat(ii).cutint, ...
                'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num)]) % post aligned by treatment time

%             concut = dat(ii).cuttime<=dat(ii).fitDat.postFitVals(1) & dat(ii).fitDat.postFitVals(1)>=dat(ii).param.t_cpd; % only pre-constricted post-treatment
%             t_concut = dat(ii).cuttime(concut);
%             i_concut = dat(ii).cutint(concut);
%             plot(t_concut-dat(ii).param.t_cpd, i_concut, ...
%                 'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num)]) % post aligned by treatment time, cut after constriction start
        end
    end
end

%% Plot axial widths vs. time
if plot_width_axial
    figure('FileName', [path '/' today '_width_axial.fig'])
    hold on
    box on
    
    for ii = 1:length(dat)
        if ~isempty(dat(ii).fitDat.preFitVals)
            plot(dat(ii).cuttime_ax - dat(ii).fitDat.preFitVals(1), dat(ii).cutdiams_ax, ...
                'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num)])
        end
    end
end

%% Plot fitted d0 value
if plot_d0
    figure('FileName', [path '/' today '_d0.fig'])
    hold on
    box on
    
    d0 = [];
    for ii = 1:length(dat)
        if ~isempty(dat(ii).fitDat) && ~isempty(dat(ii).fitDat.postRsq) && dat(ii).fitDat.postRsq > Rsq_thresh
            d0(ii) = dat(ii).fitDat.postFitVals(3);
        end
%         if ~isempty(dat(ii).fitDat) && ~isempty(dat(ii).fitDat.preRsq) && dat(ii).fitDat.preRsq > Rsq_thresh
%             d0(ii) = dat(ii).fitDat.preFitVals(3);
%         end
    end
    
    d0(d0==0) = [];
    av_d0 = mean(d0)
    hist(d0)
    text(av_d0, 3, ['\mu = ' num2str(av_d0,'%0.2f') ' \pm ' num2str(std(d0),'%0.2f') ' nm'])
end
