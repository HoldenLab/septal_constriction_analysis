% Author: Kevin Whitley
% Date created: 191106

% This script is used for fitting and plotting various data from
% constriction experiments. Input is a structure variable called alldat2.

% Version 2:
% - removed 'rawDat' stuff. Data now starts totally unfiltered (with
% 'datUnfilt', and is filtered in this script if wanted.
% - moved end-of-trace filter (to remove flat parts post-constriction) to
% filter section (was in fitting section).

% This version was used for revision 2 of the Whitley et al. 2020 FtsZ
% treadmilling paper. It's a pretty messy version, but for reproducibility
% is the version used for that submission. Future versions will be nicer.

ud.filter_dat = 1; % filter data based on fitting error, etc.
    ud.err_thresh_thick_R2 = 0.90; % thickness error threshold. error is R^2
    ud.err_thresh_diam_R2 = 0.90; % diameter error threshold. error is R^2 (new format)
    ud.err_thresh_diam = 1.1; % diameter error threshold. error is residual norm (old format)
ud.fit_dat = 1; % process and fit data to selected model for constriction
    ud.model = 'parabolic'; % constriction model to use
    ud.err_thresh_constriction_R2 = 0.8; % R^2 threshold for data to be plotted. also threshold for fitted values to be included in teff violin plots.
ud.segment_thickness = 1; % segment thickness traces
    ud.state_change = 1; % use state change to segment thicknesses.

plot_unfitted_dat = 0; % plot data traces - filtered or raw
plot_fitted_dat = 1; % plot processed traces with fitted model
    ud.alignment = 't0'; % align fitted traces. 't0' = constriction start time, 'tcpd' = compound arrival time
    plot_model = 0; % plot fitted constriction model on top of traces
    
plot_teff = 0; % make violin plots of effective constriction time (teff)
    ud.d0 = 1100; % fix d0 for calculating teff. updated 200625.
plot_perturbed = 0; % plot only traces where data was recorded DURING drug treatment
plot_intensity = 0; % plot septal intensities over time
plot_scatter = 1; % plot thickness vs. diameter scatter with division stage classification

%% Load variable dat, which has all raw (unfiltered) data from alldat2
dat=[];
for ii = 1:length(alldat2)
    for jj = 1:length(alldat2(ii).datUnfilt)
        
        if str2double(alldat2(ii).param.analysis_date) < 200321 % old format. date is mostly a guess.
            loaddat =           alldat2(ii).rawDat(jj);
            loaddat.FWHM_ax =   alldat2(ii).rawDat(jj).cutdiams_ax;
            loaddat.diam =      alldat2(ii).rawDat(jj).cutdiams;
            loaddat.time =      loaddat.cuttime;
            loaddat.time_ax =   loaddat.cuttime_ax;
        else
            loaddat =           alldat2(ii).datUnfilt(jj);
        end
        
        if isfield(alldat2,'rawDat')
            loaddat.num =       alldat2(ii).rawDat(jj).num;
        else
            loaddat.num =       alldat2(ii).datUnfilt(jj).num;
        end
        loaddat.param =         alldat2(ii).param;
        loaddat.im_date =       alldat2(ii).param.im_file(1:6);
        posind =                strfind(alldat2(ii).param.im_file, 'Pos');
        loaddat.pos =           alldat2(ii).param.im_file(posind:posind+3);
        
        dat = [dat; loaddat];
    end
end

%% Filter data frame-by-frame
if ud.filter_dat
    
    % threshold for cutting on fitted line centers. old data did not record
    % these values. using 2 sigma from mean as threshold.
    if isfield(dat(1),'linecents')
        allcents = cat(1,dat(:).linecents);
        ud.err_thresh_linecent_min = nanmean(allcents) - 2*nanstd(allcents);
        ud.err_thresh_linecent_max = nanmean(allcents) + 2*nanstd(allcents);
    end
    
    % filter each trajectory separately. would love to use a matrix or
    % table for this, but can't easily do all filters that way.
    for ii = 1:length(dat)
        
        time_dat =      dat(ii).time(:);
        thick_dat =     dat(ii).FWHM_ax(:);
        diam_dat =      dat(ii).diam(:);
        int_dat =       dat(ii).intensity(:);
        
        if str2double(dat(ii).param.analysis_date) > 200321 % new format. data starts unfiltered.
            
            thick_dat_err =     dat(ii).se_ax(:);
            diam_dat_err =      dat(ii).fiterrs(:);

            % cut on fitted center of line profile (if center close to edge
            % of image it often gives huge diameters with artificially
            % small errors).
            if isfield(dat(ii),'linecents')
                centdat = dat(ii).linecents(:);
                is_cent_ok = centdat>ud.err_thresh_linecent_min & centdat<ud.err_thresh_linecent_max;
                
                time_dat =          time_dat(is_cent_ok);
                thick_dat_err =     thick_dat_err(is_cent_ok);
                thick_dat =         thick_dat(is_cent_ok);
                diam_dat =          diam_dat(is_cent_ok);
                diam_dat_err =      diam_dat_err(is_cent_ok);
                int_dat =           int_dat(is_cent_ok);
            end
            
            % cut axial and radial times separately
            tdat_rad = time_dat;
            tdat_ax = time_dat;
            
            
            raddat2 = diam_dat;
            axdat2 = thick_dat;
            tdat2 = time_dat;
            
            
            % cut on error in thickness and radius fit
            if str2double(dat(ii).param.analysis_date) > 201018 % newer format. error using R^2
                tdat_ax = tdat_ax(thick_dat_err > ud.err_thresh_thick_R2); % R^2
                thick_dat = thick_dat(thick_dat_err > ud.err_thresh_thick_R2); % R^2
                
                tdat_rad = tdat_rad(diam_dat_err > ud.err_thresh_diam_R2); % R^2
                diam_dat = diam_dat(diam_dat_err > ud.err_thresh_diam_R2); % R^2
                int_dat = int_dat(diam_dat_err > ud.err_thresh_diam_R2);
            else % older format. error using residual squared of norms?
                tdat_rad = tdat_rad(diam_dat_err < ud.err_thresh_diam); % resnorm
                diam_dat = diam_dat(diam_dat_err < ud.err_thresh_diam); % resnorm
                int_dat = int_dat(diam_dat_err < ud.err_thresh_diam);
            end
        else % old format. data already mostly filtered.
            tdat_ax = dat(ii).time_ax(:);
            tdat_rad = time_dat;
        end
        
        % cut flat ends of traces (lingering intensity, not really
        % constriction) using first derivative
        d_s = smooth(diam_dat,10); % smooth factor was 10. now 5 (200515).
        t_s = smooth(tdat_rad,10);
        
        dDdt = (d_s(2:end)-d_s(1:end-1)) ./ (t_s(2:end)-t_s(1:end-1)); % first derivative of diameters
        dDdt_s = smooth(dDdt,3); % smoothed first derivative. factor was 10. now 3 (200515).
        dDdt_thresh = 5; % threshold for first derivative
        quartind = floor(length(dDdt_s)/4) + 1; % quarterway mark of trace (can't remove any flat parts at beginning). +1 needed to handle super short traces.
        indx = min([find(dDdt_s(quartind:end)>=dDdt_thresh, 1, 'first')+quartind-1 length(tdat_rad)]);

        tdat_rad_c = tdat_rad(1:indx); % cut times (flat ends removed)
        raddat_c = diam_dat(1:indx); % cut diameters
        ints_c = int_dat(1:indx); % cut intensities
        if ~isempty(tdat_rad_c)
            tdat_ax_c = tdat_ax(tdat_ax<tdat_rad_c(end));
            axdat_c = thick_dat(tdat_ax<tdat_rad_c(end));
        else
            tdat_ax_c = tdat_ax;
            axdat_c = thick_dat;
        end
          
        
        
        % cut flat ends of traces (lingering intensity, not really
        % constriction) using first derivative
        d_s = smooth(raddat2,10); % smooth factor was 10. now 5 (200515).
        t_s = smooth(tdat2,10);
        dDdt = (d_s(2:end)-d_s(1:end-1)) ./ (t_s(2:end)-t_s(1:end-1)); % first derivative of diameters
        dDdt_s = smooth(dDdt,3); % smoothed first derivative. factor was 10. now 3 (200515).
        
%         dDdt_s = smooth(dDdt,10);
%         halfind = floor(length(dDdt_s)/2) + 1; % +1 is extra test
%         indx = min([find(dDdt_s(halfind:end)>=dDdt_thresh, 1, 'first')+halfind length(tdat_rad)]);
        
        quartind = floor(length(dDdt_s)/4) + 1; % quarterway mark of trace (can't remove any flat parts at beginning). +1 needed to handle super short traces.
        indx = min([find(dDdt_s(quartind:end)>=dDdt_thresh, 1, 'first')+quartind-1 length(tdat2)]);
        
        tdat_rad_c2 = tdat2(1:indx); % cut times (flat ends removed)
        raddat_c2 = raddat2(1:indx); % cut diameters
        if ~isempty(tdat_rad_c2)
            tdat_ax_c2 = tdat2(tdat2<tdat_rad_c2(end));
            axdat_c2 = axdat2(tdat2<tdat_rad_c2(end));
        else
            tdat_ax_c2 = tdat2;
            axdat_c2 = axdat2;
        end
        
        
        
        
        dat(ii).cuttime_ax = tdat_ax_c;
        dat(ii).cutdiams_ax = axdat_c;
        dat(ii).cuttime_all = tdat_ax_c2;
        dat(ii).cutdiams_ax2 = axdat_c2;
        dat(ii).cutdiams2 = raddat_c2;
        dat(ii).cuttime = tdat_rad_c;
        dat(ii).cutdiams = raddat_c;
%         dat(ii).cutint = ints_c;
    end
else
    for ii = 1:length(dat)
        dat(ii).cuttime_ax = dat(ii).time(~isnan(dat(ii).FWHM_ax));
        dat(ii).cutdiams_ax = dat(ii).FWHM_ax(~isnan(dat(ii).FWHM_ax));
        dat(ii).cuttime = dat(ii).time(~isnan(dat(ii).diam));
        dat(ii).cutdiams = dat(ii).diam(~isnan(dat(ii).diam));
        dat(ii).cutint = dat(ii).intensity(~isnan(dat(ii).intensity));
    end
end

%% Plot filtered data (no fits)
if plot_unfitted_dat
    
    figure('FileName', [path '/' datestr(now,'yymmdd') '_filtered_dat.fig'])
    
    for ii = 1:length(dat)
        ax1 = subplot(211);
        hold(ax1, 'on')
        plot(ax1, dat(ii).cuttime, dat(ii).cutdiams, 'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num)])
        ylim([200 1200])
        
        ax2 = subplot(212);
        hold(ax2,'on')
        plot(dat(ii).cuttime_ax, dat(ii).cutdiams_ax, 'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num)])
        
        % partition thickness traces, plot separately
        %         cond = detectZRcondensation(dat(ii).cutdiams_ax, 'MinSpatialDist', 50);
        %         plot(ax2,dat(ii).cuttime_ax(cond==1), dat(ii).cutdiams_ax(cond==1), 'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num)], 'Color', 'k')
        %         plot(ax2,dat(ii).cuttime_ax(cond==2), dat(ii).cutdiams_ax(cond==2), 'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num)], 'Color', 'r')
        ylim([200 1000])
    end
    linkaxes([ax1 ax2],'x')
    
    tv = 0:1000;
    plot(ones(length(tv),1)*dat(ii).param.t_cpd, tv, 'r')
    
    title(ax1, 'Filtered data')
    ylabel(ax1, 'Diameter (nm)')
    xlabel(ax2, 'Time (min)')
    ylabel(ax2, 'Thickness (nm)')
    set(gcf, 'UserData', ud)
end

%% Process and fit data to selected model for constriction
if ud.fit_dat
    
    if plot_fitted_dat
        figure('FileName', [path '/' datestr(now,'yymmdd') '_fit_dat.fig'])
        ax1 = subplot(211);
        hold(ax1, 'on')
        box(ax1, 'on')
        ax2 = subplot(212);
        hold(ax2, 'on')
        box(ax2, 'on')
%         ax3 = subplot(313);
%         hold(ax3, 'on')
%         box(ax3, 'on')
    end
    
    diffp=[]; diffl=[]; dax_precons=[]; dax_postcons=[]; postcon_num=0; precons_decond=[]; precons_cond=[]; std_precons=[];
    for ii = 1:length(dat)
        
        % 1. prepare data for fitting
        
        % make sure these are both column vectors
        tfit = dat(ii).cuttime(:);
        dfit = dat(ii).cutdiams(:);
        tax_c = dat(ii).cuttime_ax(:);
        dax_c = dat(ii).cutdiams_ax(:);

        % Split data into pre- and post-treatment
        % exception: if param.t_cpd=0, assume all data is untreated (or
        % 'pre-treated')
        if dat(ii).param.t_cpd == 0
            dat(ii).preTimeCut = tfit;
            dat(ii).preDiamCut = dfit;
            dat(ii).postTimeCut = [];
            dat(ii).postDiamCut = [];
        else
            dat(ii).preTimeCut = tfit(tfit<dat(ii).param.t_cpd);
            dat(ii).preDiamCut = dfit(tfit<dat(ii).param.t_cpd);
            dat(ii).postTimeCut = tfit(tfit>=dat(ii).param.t_cpd);
            dat(ii).postDiamCut = dfit(tfit>=dat(ii).param.t_cpd);
        end
        
        % 2. fit the data
        
        % fit pre-treatment data
        if length(dat(ii).preTimeCut) >= 5 % needs to have enough points for decent fit
            [cfit_pre, rn_pre, res, fcnp] = fit_constriction_data(dat(ii).preTimeCut, dat(ii).preDiamCut, ud.model); % parabolic fit
%             [cfit_prel, rn_prel, resl, fcnl] = fit_constriction_data(dat(ii).preTimeCut, dat(ii).preDiamCut, 'linear'); % linear fit
            
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
        
        dat(ii).fcn = fcnp;
        dat(ii).fitDat.num = ii;
        
        dat(ii).fitDat.preFitVals = cfit_pre;
        dat(ii).fitDat.preRsq = Rsq_pre;
        dat(ii).fitDat.preRsq_con = Rsq_pre_con;
        dat(ii).fitDat.n_pre = n_pre;
        dat(ii).fitDat.XSq = Xsq_pre;
        dat(ii).fitDat.XSq_pre_con = Xsq_pre_con;

        % fit post-treatment data
        if length(dat(ii).postTimeCut) >= 5 && any(dat(ii).postTimeCut <= dat(ii).param.t_cpd+1) % bit shoddy, but works
            [cfit_post, rn_post, res, fcnp] = fit_constriction_data(dat(ii).postTimeCut, dat(ii).postDiamCut, ud.model);
            
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

        dat(ii).fitDat.postFitVals = cfit_post;
        dat(ii).fitDat.postRsq = Rsq_post;
        dat(ii).fitDat.postRsq_con = Rsq_post_con;
        dat(ii).fitDat.n_post = n_post;

        % partition thickness traces into decondensed and condensed
        if ud.segment_thickness
            
            seg_t={}; seg_d={}; clr={}; dur=[];
%             if Rsq_pre_con > Rsq_thresh && n_pre>5
            if n_pre>5 % need Rsq? 210120

                % using state change
                if ud.state_change
%                     [inds, switchPt, muse] = detectZRcondensation(dax_c, 'MinSpatialDist', 50, 'Statistic', 'std'); % from point change
                    
                    
                    % state change analysis only on pre-constriction data
                    % 210121
                    [inds, switchPt, muse] = detectZRcondensation(dax_c(tax_c<cfit_pre(1)), 'MinSpatialDist', 50, 'Statistic', 'std'); % from point change
                    
                    
                 
                    inds(inds==2) = 0; % for consistency with threshold case
                    
                    
%                     val_end = inds(end);
                    inds = [inds; ones(length(dax_c)-length(inds),1)]; % add all constricting points too 210202
                    
                    
                    
                    transind = [0; switchPt(switchPt~=-1); length(inds)];
                    
                    if switchPt ~= -1
                        dat(ii).fitDat.t_condense = tax_c(switchPt);
                    else
                        dat(ii).fitDat.t_condense = 0;
                    end
                else % using threshold
                    %                 thick_cut = 422; % [nm] mean + 3*sigma for 'condensed' state of wt FtsZ-GFP
                    %                 thick_cut = 389;% [nm] mean + 2*sigma for 'condensed' state of wt FtsZ-GFP
                    thick_cut = 356;% [nm] mean + 1*sigma for 'condensed' state of wt FtsZ-GFP
                    %                 dax_c_sm = smooth(dax_c, 3);
                    dax_c_sm = dax_c;
                    inds = dax_c_sm < thick_cut;
                    trans = inds(1:end-1) - inds(2:end);
                    transind = [0; find(trans==-1); find(trans==1); length(inds)];
                    transind = sort(transind);
                end
                
                % segmentation
                for kk = 1:length(transind)-1
                    seg_t{kk} = tax_c(transind(kk)+1:transind(kk+1));
                    seg_d{kk} = dax_c(transind(kk)+1:transind(kk+1));
                    if ~isempty(tax_c) && transind(kk)+1<tax_c(end) && ~isempty(inds) && inds(transind(kk)+1)==0 % decondensed
                        clr{kk} = 'r';
                        
%                         precons_decond = [precons_decond; seg_d{kk}];
                    elseif ~isempty(tax_c) % condensed
                        clr{kk} = 'k';
                        
                        % find every consecutive stretch of 'condensed'
                        % state prior to constriction, measure length of
                        % segment
                        precons = seg_t{kk}-cfit_pre(1);
                        precons = precons(precons<=0);
                        dur(kk) = length(precons);
                        %  dur(kk) = length(seg_t{kk});
%                         precons_cond = [precons_cond; seg_d{kk}];
                    end
                end
                
                dat(ii).cond_dur = max([0 dur]);
                
                allprecons = dat(ii).preTimeCut-cfit_pre(1);
                allprecons = allprecons(allprecons<=0); % need to be sure there are actually points there to compare
                if length(allprecons)>=1 && ~isempty(dat(ii).cond_dur) && dat(ii).cond_dur>=1 % condensation for 1 point pre-constriction
                    dat(ii).didcond = 1;
                elseif length(allprecons)>=1 && ~isempty(dat(ii).cond_dur) && dat(ii).cond_dur<1 % no condensation for 1 point pre-constriction (but there was 1 data point)
                    dat(ii).didcond = 0;
                end
                
                postcon_num = postcon_num+1;
                dax_precons = [dax_precons; dax_c(tax_c<cfit_pre(1))];
                dax_postcons = [dax_postcons; dax_c(tax_c>=cfit_pre(1))];
            elseif n_pre>5 && ~isempty(tax_c)
                transind = [0 tax_c(end)];
                seg_t{1} = tax_c;
                seg_d{1} = dax_c;
                clr{1} = 'b';
            else
                transind = 0;
            end
        end
        
        % 3. plot the data with fits
        
        if plot_fitted_dat
            
            % plot data pre-treatment
            if Rsq_pre_con > ud.err_thresh_constriction_R2 && n_pre>5
                
                if strcmp(ud.alignment, 'tcpd')
                    shift = dat(ii).param.t_cpd; % define t=0 as treatment time
                elseif strcmp(ud.alignment, 't0')
                    shift = cfit_pre(1); % define t=0 as constriction start time
                end
                
                % hold(ax1,'off')
                % hold(ax2,'off')
                
                % plot diameters, shifted by either treatment time or
                % constriction start time
                plot(ax1, dat(ii).preTimeCut-shift, dat(ii).preDiamCut, 'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num) '\_pre'])
                
                if plot_model
                    tvf = dat(ii).preTimeCut(1):0.1:dat(ii).preTimeCut(end);
                    plot(ax1, tvf-shift, fcnp(cfit_pre,tvf), 'k', 'DisplayName', '')
                end
                
                %                 circ = dat(ii).cutdiams*pi; % circumferences
                %                 plot(ax2, dat(ii).cuttime-shift, dat(ii).cutint, 'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num) '\_pre']) % density
                
%                 plot(ax2, tax_c-shift, dat(ii).cutstd_ax, 'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num) '\_pre'])
                
%                 ind_tcon = find(abs(tax_c-cfit_pre(1)) < 0.5, 1, 'first');
%                 if ~isempty(ind_tcon)
%                     dense_ax = dat(ii).cutint_ax ./ dax_c;
%                     ax_int_con = dense_ax(ind_tcon);
%                     dense_ax = dense_ax / ax_int_con;
%                     plot(ax2, tax_c-shift, dense_ax, 'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num) '\_pre']) % density
%                 end
                
%                 wind = 7;
%                 [~, ~, pval] = ttesttrace2(tax_c, dax_c, 0.3, wind);
%                 plot(ax2, tax_c(1:end-wind)-shift, -log10(pval), 'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num) '\_pre'])

                if ud.segment_thickness
                    % plot segmented thickness traces
                    for kk = 1:length(transind)-1
                        if ~isempty(tax_c)
                            plot(ax2, seg_t{kk}-shift, seg_d{kk}, 'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num) '\_pre'], 'Color', clr{kk})
                            %                     hold(ax2,'on')
                        end
                    end
                %                 plot(ax2, tax_c(cond==1)-shift, dax_c(cond==1), 'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num) '\_pre'], 'linew', 1, 'Color', 'k')
                %                 plot(ax2, tax_cond-shift, dax_cond, 'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num) '\_pre'], 'linew', 1, 'Color', 'k')
                %                 hold(ax2,'on')
                %                 plot(ax2, tax_c(cond==2)-shift, dax_c(cond==2), 'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num) '\_pre'], 'linew', 1, 'Color', 'r')
                %                 plot(ax2, tax_decond-shift, dax_decond, 'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num) '\_pre'], 'linew', 1, 'Color', 'r')
                else
%                     dax_precons = [dax_precons; dax_c(tax_c<cfit_pre(1))];
                    dax_precons = dax_c(tax_c<cfit_pre(1));
                    std_precons(ii) = std(dax_precons);
                    plot(ax2, tax_c-shift, dax_c, 'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num) '\_pre'])
                end

                % calculate difference between 2B disappearance time and
                % predicted end of constriction. Since d0 is not constant
                % for each trace, can't use simple formula for end time of
                % constriction.
%                 gone2b = dat(ii).cuttime(end);
%                 tvf2 = dat(ii).preTime(1):0.01:dat(ii).preTime(end)+50;
%                 pararray = fcnp(cfit_pre, tvf2);
%                 parzero = find(~pararray, 1, 'first');
%                 teffp = tvf2(parzero); % can't use fzero because lot of zeros
%                 teffl = fzero(@(x)fcnl(cfit_prel,x), cfit_pre(1));
                
%                 diffp = [diffp; teffp-gone2b];
%                 diffl = [diffl; teffl-gone2b];
            end
            
            % plot data post-treatment (if exists)
            if Rsq_post_con > ud.err_thresh_constriction_R2 && n_post>5
%                 plot(ax1, dat(ii).postTimeCut-shift, dat(ii).postDiamCut, 'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num) '\_post'], 'linew', 1)
                
                if plot_model
                    tvf = dat(ii).postTimeCut(1):0.1:dat(ii).postTimeCut(end);
                    plot(ax1, tvf-shift, fcnp(cfit_post,tvf), 'r', 'DisplayName', '')
                end
            end
            
        end
    end
    linkaxes([ax1 ax2], 'x')
    set(gcf, 'UserDat', ud)
end

%% Make violin plots (with DABEST if possible) of effective constriction times
if plot_teff

    teff_pre = []; teff_post = []; cat_pre = {}; cat_post = {};
    for ii = 1:length(dat)

        if ~isempty(dat(ii).fitDat) && dat(ii).fitDat.preRsq_con > ud.err_thresh_constriction_R2 && dat(ii).fitDat.n_pre>5
            teff_pre = [teff_pre ud.d0^2 ./ dat(ii).fitDat.preFitVals(2)];
            cat_pre = [cat_pre 'Untreated'];
%         else
%             teff_pre(ii) = [];
        end
        
        if ~isempty(dat(ii).fitDat) && dat(ii).fitDat.postRsq_con > ud.err_thresh_constriction_R2 && dat(ii).fitDat.n_post>5
            teff_post = [teff_post ud.d0^2 ./ dat(ii).fitDat.postFitVals(2)];
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
        figure('FileName', [path '/' datestr(now,'yymmdd') '_teff.fig'])
        
        violinplot(teff_pre)
        
        ylabel('t_{eff} (min)')
    end
    set(gcf, 'UserDat', ud)
end

%% Plot only traces where data was recorded DURING drug treatment
if plot_perturbed
    figure('FileName', [path '/' datestr(now,'yymmdd') '_perturbed.fig'])
    hold on
    box on
    
    for ii = 1:length(dat)
        if any(dat(ii).cuttime<dat(ii).param.t_cpd) && any(dat(ii).cuttime>dat(ii).param.t_cpd)
            plot(dat(ii).cuttime, dat(ii).cutdiams, 'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num)])
        end
    end
end

%% Plot septal intensites vs. time
if plot_intensity
    figure('FileName', [path '/' datestr(now,'yymmdd') '_pre.fig'])
    subplot(212)
    h1 = gca;
    hold on
    subplot(211)
    g1 = gca;
    hold on
    
    figure('FileName', [path '/' datestr(now,'yymmdd') '_post.fig'])
    subplot(212)
    h2 = gca;
    hold on
    subplot(211)
    g2 = gca;
    hold on
    
    for ii = 1:length(dat)
        
        if ~any(dat(ii).cuttime>dat(ii).param.t_cpd) && ~isempty(dat(ii).fitDat.preFitVals) && dat(ii).fitDat.preRsq_con > ud.err_thresh_constriction_R2 && dat(ii).fitDat.n_pre>5 % only before
%         if any(dat(ii).cuttime<dat(ii).param.t_cpd) && any(dat(ii).cuttime>dat(ii).param.t_cpd) % only perturbed
%         if isfield(dat(ii).fitDat,'preFitVals') && ~isempty(dat(ii).fitDat.preFitVals)

            plot(g1, dat(ii).cuttime-dat(ii).fitDat.preFitVals(1), dat(ii).cutdiams, ...
                'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num)]) % post aligned by treatment time

            plot(h1, dat(ii).cuttime-dat(ii).fitDat.preFitVals(1), dat(ii).cutint, ...
                'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num)]) % pre aligned by constriction time
        end
        
%         if any(dat(ii).cuttime>dat(ii).param.t_cpd) && ~isempty(dat(ii).fitDat.postFitVals) && dat(ii).fitDat.postRsq_con > Rsq_thresh && dat(ii).fitDat.n_post>5 % only after
%         if any(dat(ii).cuttime>dat(ii).param.t_cpd) && any(dat(ii).cuttime>dat(ii).param.t_cpd) && isempty(dat(ii).fitDat.postFitVals) && ~any(dat(ii).cutdiams<700)% only after, raw data
%         if any(dat(ii).cuttime==dat(ii).param.t_cpd) && dat(ii).cutdiams(dat(ii).cuttime==dat(ii).param.t_cpd)>900 % only after, mature rings only
        if dat(ii).cutdiams(dat(ii).cuttime==dat(ii).param.t_cpd)>900 % only after, mature rings only
            
%             plot(dat(ii).cuttime-dat(ii).fitDat.postFitVals(1), dat(ii).cutint, ...
%                 'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num)]) % post aligned by constriction time

            plot(g2, dat(ii).cuttime-dat(ii).fitDat.preFitVals(1), dat(ii).cutdiams, ...
                'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num)]) % post aligned by treatment time
            
            circ = dat(ii).cutdiams*pi; % circumferences
            
            plot(h2, dat(ii).cuttime-dat(ii).fitDat.preFitVals(1), dat(ii).cutint./circ, ...
                'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num)]) % post aligned by treatment time

%             concut = dat(ii).cuttime<=dat(ii).fitDat.postFitVals(1) & dat(ii).fitDat.postFitVals(1)>=dat(ii).param.t_cpd; % only pre-constricted post-treatment
%             t_concut = dat(ii).cuttime(concut);
%             i_concut = dat(ii).cutint(concut);
%             plot(t_concut-dat(ii).param.t_cpd, i_concut, ...
%                 'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num)]) % post aligned by treatment time, cut after constriction start
        end
    end
end

if plot_scatter
    
    figure('FileName', [path '/' datestr(now,'yymmdd') '_scatter.fig'])
    hold on
    box on
    
    cons = []; nasc = []; mats = []; stds_ax = []; twostate = []; stds_mat = []; stds_nasc = [];
    for ii = 1:length(dat)
        if ~isempty(dat(ii).fitDat.preFitVals) && dat(ii).fitDat.preRsq_con > ud.err_thresh_constriction_R2 && dat(ii).fitDat.preFitVals(1)>dat(ii).cuttime_all(1)-5

            if ud.segment_thickness
                if isfield(dat(ii).fitDat,'t_condense') && ~isempty(dat(ii).fitDat.t_condense)
                    inds_nascent = dat(ii).cuttime_all < dat(ii).fitDat.t_condense;
                    inds_mature = dat(ii).cuttime_all >= dat(ii).fitDat.t_condense & dat(ii).cuttime_all < dat(ii).fitDat.preFitVals(1);
                else
                    inds_nascent = false(length(dat(ii).cuttime_all),1);
                    inds_mature = dat(ii).cuttime_all < dat(ii).fitDat.preFitVals(1);
                end
                
                if any(inds_nascent) && any(inds_mature)
                    twostate = [twostate; 1]; % showed two-state behavior
                else
                    twostate = [twostate; 0]; % did not show two-state behavior
                end
                
                nasc = [nasc; dat(ii).cutdiams_ax2(inds_nascent) dat(ii).cutdiams2(inds_nascent)];
                mats = [mats; dat(ii).cutdiams_ax2(inds_mature) dat(ii).cutdiams2(inds_mature)];
                cons = [cons; dat(ii).cutdiams_ax2(~(inds_nascent | inds_mature)) dat(ii).cutdiams2(~(inds_nascent | inds_mature))];
                
                stds_ax = [stds_ax; std([dat(ii).cutdiams_ax2(inds_nascent); dat(ii).cutdiams_ax2(inds_mature)])];
                stds_nasc = [stds_nasc; std([dat(ii).cutdiams_ax2(inds_nascent)])];
                stds_mat = [stds_mat; std([dat(ii).cutdiams_ax2(inds_mature)])];
            else
                inds_precon = dat(ii).cuttime_all < dat(ii).fitDat.preFitVals(1);
                
                mats = [mats; dat(ii).cutdiams_ax2(inds_precon) dat(ii).cutdiams2(inds_precon)];
                cons = [cons; dat(ii).cutdiams_ax2(~inds_precon) dat(ii).cutdiams2(~inds_precon)];
            end
        end
    end
    
    plot(cons(:,2), cons(:,1), 'ob')
    plot(mats(:,2), mats(:,1), 'or')
    if ud.segment_thickness
        plot(nasc(:,2), nasc(:,1), 'og')
    end
    xlim([200 1000])
    ylim([200 1000])
    set(gcf, 'UserDat', ud)
end