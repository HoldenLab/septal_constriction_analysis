% Author: Kevin Whitley
% Date created: 191106

% This script is used for fitting and plotting various data from
% constriction experiments. Input is a structure variable called alldat2,
% which must contain the field rawDat.

fit_dat = 1; % process and fit data to selected model for constriction
model = 'parabolic'; % constriction model to use
Rsq_thresh = 0.8; % R^2 threshold for data to be plotted. also threshold for fitted values to be included in teff violin plots.
d0 = 1100; % fix d0 for calculating teff. updated 200625.

plot_raw_dat = 0; % plot raw data traces
plot_fit_dat = 0; % plot processed traces with fitted model
plot_teff = 1; % make violin plots of effective constriction time (teff)
plot_d_vs_alph = 0; % plot diameter at time of drug treatment vs. constriction rate
plot_perturbed = 0; % plot only traces where data was recorded DURING drug treatment
plot_intensity = 0; % plot septal intensities over time
plot_width_axial = 0; % plot axial widths (thickness) over time
plot_d0 = 0; % plot fitted d0 (pretty useless)

%% Load variable dat, which has all raw data from alldat2
dat = [];
for ii = 1:length(alldat2)
    for jj = 1:length(alldat2(ii).rawDat)
        rawdat = alldat2(ii).rawDat(jj);
        
        rawdat.param = alldat2(ii).param;
        rawdat.im_date = alldat2(ii).param.im_file(1:6);
        
        posind = strfind(alldat2(ii).param.im_file, 'Pos');
        rawdat.pos = alldat2(ii).param.im_file(posind:posind+3);
        
        dat = [dat; rawdat];
    end
end

%% Plot all raw data traces
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

%% Process and fit data to selected model for constriction
if fit_dat
    
    if plot_fit_dat
        figure('FileName', [path '/' today '_fit_dat.fig'])
        hold on
        box on
    end
    
    diffp=[]; diffl=[];
    for ii = 1:length(dat)
        
        % 1. prepare data for fitting
        
        % make sure these are both column vectors
        tfit = dat(ii).cuttime(:);
        dfit = dat(ii).cutdiams(:);

        % trim to remove flat ends of traces based on first derivative
        % (this is not constriction).
        d_s = smooth(dfit,5); % smooth factor was 10. now 5 (200515).
        t_s = smooth(tfit,5);
        dDdt = (d_s(2:end)-d_s(1:end-1)) ./ (t_s(2:end)-t_s(1:end-1)); % first derivative of diameters
        dDdt_s = smooth(dDdt,3); % smoothed first derivative. factor was 10. now 3 (200515).
        dDdt_thresh = 0; % threshold for first derivative
        halfind = floor(length(dDdt_s)/2) + 1; % halfway mark of trace (can't remove any flat parts at beginning). +1 needed to handle super short traces.
        indx = min([find(dDdt_s(halfind:end)>=dDdt_thresh, 1, 'first')+halfind-1 length(tfit)]);
        
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
        
        % Split data into pre- and post-treatment
        % exception: if param.t_cpd=0, assume all data is untreated (or
        % 'pre-treated')
        if dat(ii).param.t_cpd == 0
            dat(ii).preTime = tfit;
            dat(ii).preDiam = dfit;
            dat(ii).postTime = [];
            dat(ii).postDiam = [];
            
            dat(ii).preTimeCut = tfit_c;
            dat(ii).preDiamCut = dfit_c;
            dat(ii).postTimeCut = [];
            dat(ii).postDiamCut = [];
        else
            dat(ii).preTime = tfit(tfit<dat(ii).param.t_cpd);
            dat(ii).preDiam = dfit(tfit<dat(ii).param.t_cpd);
            dat(ii).postTime = tfit(tfit>=dat(ii).param.t_cpd);
            dat(ii).postDiam = dfit(tfit>=dat(ii).param.t_cpd);
            
            dat(ii).preTimeCut = tfit_c(tfit_c<dat(ii).param.t_cpd);
            dat(ii).preDiamCut = dfit_c(tfit_c<dat(ii).param.t_cpd);
            dat(ii).postTimeCut = tfit_c(tfit_c>=dat(ii).param.t_cpd);
            dat(ii).postDiamCut = dfit_c(tfit_c>=dat(ii).param.t_cpd);
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
            fcn = [];
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
        
        if plot_fit_dat
            if Rsq_pre_con > Rsq_thresh && n_pre>5
%                 plot(dat(ii).preTime-dat(ii).param.t_cpd, dat(ii).preDiam, 'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num) '\_pre'], 'linew', 1)
                plot(dat(ii).preTimeCut-dat(ii).param.t_cpd, dat(ii).preDiamCut, 'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num) '\_pre'], 'linew', 1)
%                 plot(dat(ii).preTime-cfit_pre(1), dat(ii).preDiam, 'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num) '\_pre'], 'linew', 1)

                tvf = dat(ii).preTime(1):0.1:dat(ii).preTime(end);
                plot(tvf-dat(ii).param.t_cpd, fcnp(cfit_pre,tvf), 'k', 'DisplayName', '')
%                 plot(tvf-dat(ii).param.t_cpd, fcnl(cfit_prel,tvf), 'k', 'DisplayName', '')
%                 plot(tvf-cfit_pre(1), fcn(cfit_pre,tvf), 'k', 'DisplayName', '')

                % calculate difference between 2B disappearance time and
                % predicted end of constriction. Since d0 is not constant
                % for each trace, can't use simple formula for end time of
                % constriction.
                gone2b = dat(ii).cuttime(end);
                tvf2 = dat(ii).preTime(1):0.01:dat(ii).preTime(end)+50;
                pararray = fcnp(cfit_pre, tvf2);
%                 teffp = d0^2 ./ cfit_pre(2);
                parzero = find(~pararray, 1, 'first');
                teffp = tvf2(parzero); % can't use fzero because lot of zeros
%                 teffl = -d0 / cfit_prel(2);
                teffl = fzero(@(x)fcnl(cfit_prel,x), cfit_pre(1));
                
                diffp = [diffp; teffp-gone2b];
                diffl = [diffl; teffl-gone2b];
            end
            if Rsq_post_con > Rsq_thresh && n_post>5
%                 plot(dat(ii).postTime-dat(ii).param.t_cpd, dat(ii).postDiam, 'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num) '\_post'], 'linew', 1)
%                 plot(dat(ii).cuttime-dat(ii).param.t_cpd, dat(ii).cutdiams, 'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num) '\_post'], 'linew', 1)

                tvf = dat(ii).postTime(1):0.1:dat(ii).postTime(end);
%                 plot(tvf-dat(ii).param.t_cpd, fcnp(cfit_post,tvf), 'r', 'DisplayName', '')
                teff = d0^2 ./ cfit_post(2);
            end
        end
    end
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
        plot(dat(ii).cuttime_ax, dat(ii).cutdiams_ax, ...
            'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num)])
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
