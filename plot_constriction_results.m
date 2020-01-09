% Author: Kevin Whitley
% Date created: 191106

% This script is used for analysing and plotting the results from multiple
% FOVs of constriction experiments.

fit_dat = 1;

plot_raw_dat = 0;
plot_fit_dat = 0;
plot_filter_dat = 0;
plot_t_con = 0;
plot_alph = 0;
plot_d_vs_alph = 0;
plot_perturbed = 1;

chiSq_thresh = 6e4;
Rsq_thresh = 0.8;
% Rsq_thresh = 0;
t0_thresh = 1000; % [min] standard error from fit
alph_thresh = 6e11; % [nm^2/min] standard error from fit

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

% tcut = [0 50];
% dat = cutstruct(dat, 'cuttime', tcut);

if fit_dat
    ud = [];
    figure
    hold on
    box on
    filts1 = []; filts2 = []; filts3 = []; filts4 = []; rn = []; se = []; Rsq = [];
    for ii = 1:length(dat)
        
        % cut anything starting after treatment time
        if isempty(dat(ii).cuttime) || dat(ii).cuttime(1) > dat(ii).param.t_cpd
            continue
        end
        
        ffit = dat(ii).cuttime'; % [min]
        if mean(dat(ii).cutdiams)<50 % earlier data was in pixels rather than nm
            lfit = dat(ii).cutdiams * 65; % [nm]
        else
            lfit = dat(ii).cutdiams;
        end

        filts1(ii) = ii;
        if length(ffit)<3 % too short - exclude
            continue
        end
        
        filts2(ii) = ii;
        if length(dat(ii).cuttime) ~= length(dat(ii).cutdiams) % not sure why this would be the case, but happened once
            continue
        end
        
        filts3(ii) = ii;
        if ~any(lfit < 400) % doesn't constrict much - fits nearly always unreliable
            continue
        end
        
        % trim to remove flat ends of traces based on first derivative
        % (this is not constriction).
        w_s = smooth(lfit,10);
        t_s = smooth(dat(ii).cuttime,10);
        
        bins = t_s(2:end)-t_s(1:end-1);
        w1d = (w_s(2:end)-w_s(1:end-1))./bins; % first derivative
        w1d_s = smooth(w1d,10);
        t1 = t_s(1:end-1) + (t_s(2:end)-t_s(1:end-1))/2;
        
        w1d_thresh = 0;
        
        halfind = floor(length(w1d_s)/2);
        indx = min([find(w1d_s(halfind:end)>=w1d_thresh, 1, 'first')+halfind length(ffit)]);
        
%         bins2 = bins(1:end-1) + (bins(2:end)-bins(1:end-1))/2;
%         w2d = (w1d_s(2:end)-w1d_s(1:end-1))./bins2; % second derivative
%         t2 = t1(1:end-1) + (t1(2:end)-t1(1:end-1))/2;
        
%         w2d_thresh = 10;
        
%         w2d_c = w2d(2:end-1);
%         t2_c = t2(2:end-1);
%         indx = min([find(w2d_c>=w2d_thresh, 1, 'first') length(ffit)]);
        
        ffit_c = ffit(1:indx);
        lfit_c = lfit(1:indx);
        
        % HACK to exclude everything pre-treatment
%         lfit_c = lfit_c(ffit_c>dat(ii).param.t_cpd);
%         ffit_c = ffit_c(ffit_c>dat(ii).param.t_cpd);
        
        filts4(ii) = ii;
        if length(ffit_c)<8 % too short - don't bother trying to fit
            continue
        end

        [cfit, rn(ii), se(ii,:), fcn] = fit_constriction_data(ffit_c, lfit_c, dat(ii).param.model);
        
        ud.fcn = fcn;
        ud.fitDat(ii).num = ii;
        ud.fitDat(ii).time = ffit;
        ud.fitDat(ii).width = lfit;
        ud.fitDat(ii).constFitVals = cfit;
        
        dat(ii).fcn = fcn;
        dat(ii).fitDat.num = ii;
        dat(ii).fitDat.time = ffit;
        dat(ii).fitDat.width = lfit;
        dat(ii).fitDat.constFitVals = cfit;
        
        sst(ii) = sum((lfit_c - mean(lfit_c)).^2);
        Rsq = 1 - rn(ii)/sst(ii);
        dat(ii).fitDat.Rsq = Rsq;
        
        chisq = rn(ii)/sqrt(length(ffit));
        ud.fitDat(ii).constFitChiSq = chisq;
        
        dat(ii).fitDat.constFitChiSq = chisq;
        
%         if chisq < chiSq_thresh && se(ii,1) < t0_thresh && se(ii,2) < alph_thresh
        if Rsq > Rsq_thresh
            plot(ffit-dat(ii).param.t_cpd, lfit, 'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num)], 'linew', 1)
%             hold on

            tvf = ffit(1):0.1:ffit(end);
            plot(tvf-dat(ii).param.t_cpd, fcn(cfit,tvf), 'k', 'DisplayName', '')
%             hold off
            
            teff = cfit(3)^2 ./ cfit(2);
        end
    end
end

if plot_raw_dat
    figure
    hold on
    box on
    
    for ii = 1:length(dat)
        plot(dat(ii).cuttime, dat(ii).cutdiams, ...
            'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num)])
    end
    tv = 0:1200;
    plot(ones(length(tv),1)*dat(ii).param.t_cpd, tv, 'r')
    
    xlabel('Time (min)')
    ylabel('d_{septum} (nm)')
    title('Raw data')
end

if plot_fit_dat
    figure
    hold on
    box on
    
%     fov = alldat.fov(jj);
    for ii = 1:length(dat)
        plot(dat.fitDat(ii).time, dat.fitDat(ii).width, 'DisplayName', num2str(dat.fitDat(ii).num))
        
        tvf = fov.fitDat(ii).time(1):0.1:fov.fitDat(ii).time(end);
        plot(tvf, fov.fcn(fov.fitDat(ii).constFitVals,tvf), 'k', 'DisplayName', '')
    end
    tv = 0:1200;
    plot(ones(length(tv),1)*fov.param.t_cpd, tv, 'r')
    
    xlabel('Time (min)')
    ylabel('d_{septum} (nm)')
    title(['Fitted, unfiltered data for ' strrep(fov.fname,'_','\_')])
end

if plot_filter_dat
    figure
    hold on
    box on
    
    for jj = 1:length(alldat.fov)
%         fov = alldat.fov(jj);
        fov = alldat2(jj);
        for ii = 1:length(fov.fitDat)
            if fov.fitDat(ii).constFitChiSq < chiSq_thresh
                plot(fov.fitDat(ii).time, fov.fitDat(ii).width, 'DisplayName', num2str(fov.fitDat(ii).num), 'linew', 1)

                tvf = fov.fitDat(ii).time(1):0.1:fov.fitDat(ii).time(end);
                plot(tvf, fov.fcn(fov.fitDat(ii).constFitVals,tvf), 'k', 'DisplayName', '')
            end
        end
    end
    tv = 0:1200;
    plot(ones(length(tv),1)*fov.param.t_cpd, tv, 'r')
    
    xlabel('Time (min)')
    ylabel('d_{septum} (nm)')
    title(['Fitted, filtered data for ' strrep(fov(1).param.tracks_file,'_','\_')])
end

if plot_t_con
    figure('FileName', [path '/' today '_t_con.fig'])
    hold on
    box on
    
    t_con = []; 
    for jj = 1:length(alldat.fov)
        fov = alldat.fov(jj);
        for ii = 1:length(fov.fitDat)
            if fov.fitDat(ii).constFitChiSq < chiSq_thresh
                t_con = [t_con; fov.fitDat(ii).time(end) - fov.fitDat(ii).constFitVals(1) fov.fitDat(ii).time(end)];
            end
        end
    end

    t_con(t_con(:,1)<1,:) = [];
    
    t_con_pre = t_con(t_con(:,2)<fov.param.t_cpd,1);
    t_con_post = t_con(t_con(:,2)>fov.param.t_cpd,1);
    tpr = histogram(t_con_pre,5:2:40);
    tpo = histogram(t_con_post,5:2:40);
    set(tpr,'facea',0.5)
    set(tpo,'facea',0.5)
    xlabel('t_{constriction} (min)')
    ylabel('Counts')
end

if plot_alph
    figure('FileName', [path '/' today '_alpha.fig'])
    hold on
    box on
    
    alph = [];
    for ii = 1:length(dat)
%         if ~isempty(dat(ii).fitDat) && dat(ii).fitDat.constFitChiSq < chiSq_thresh  && se(ii,1) < t0_thresh && se(ii,2) < alph_thresh
        if ~isempty(dat(ii).fitDat) && dat(ii).fitDat.Rsq > Rsq_thresh
            teff = dat(ii).fitDat.constFitVals(3)^2 ./ dat(ii).fitDat.constFitVals(2);
            alph = [alph; dat(ii).fitDat.constFitVals(2) dat(ii).fitDat.time(end) dat(ii).param.t_cpd teff];
        end
    end
    
    teff_pre = alph(alph(:,2)<alph(:,3),4);
    teff_post = alph(alph(:,2)>alph(:,3),4);
    alph_pre = alph(alph(:,2)<alph(:,3),1);
    alph_post = alph(alph(:,2)>alph(:,3),1);
    
    apr = histogram(alph_pre,1e4:5e3:1e5);
    apo = histogram(alph_post,1e4:5e3:1e5);
    set(apr,'facea',0.5)
    set(apo,'facea',0.5)
    xlabel('\alpha (nm^2/min)')
    ylabel('Counts')
    
    alph_pre = [alph_pre; ones(length(alph_post)-length(alph_pre),1)*NaN];
    alph_post = [alph_post; ones(length(alph_pre)-length(alph_post),1)*NaN];
    teff_pre = [teff_pre; ones(length(teff_post)-length(teff_pre),1)*NaN];
    teff_post = [teff_post; ones(length(teff_pre)-length(teff_post),1)*NaN];
    
    figure
    violinplot([teff_pre teff_post])
end

if plot_d_vs_alph
    figure('FileName', [path '/' today '_alpha.fig'])
    hold on
    box on
    
    d_cpd = []; alph_cpd = [];
    for ii = 1:length(dat)
        if ~isempty(dat(ii).fitDat) && any(dat(ii).fitDat.time<dat(ii).param.t_cpd) && any(dat(ii).fitDat.time>dat(ii).param.t_cpd) && dat(ii).fitDat.constFitChiSq < chiSq_thresh
            rng = 1; % [min]
            ind = find(dat(ii).fitDat.time>=dat(ii).param.t_cpd-rng & dat(ii).fitDat.time<=dat(ii).param.t_cpd+rng);
            d_all = dat(ii).fitDat.width(ind);
            d_cpd = [d_cpd; mean(d_all)];
            teff = dat(ii).fitDat.constFitVals(3)^2 ./ dat(ii).fitDat.constFitVals(2);
            alph_cpd = [alph_cpd; dat(ii).fitDat.constFitVals(2) teff];
        end
        
        if ~isempty(dat(ii).fitDat) && ~isempty(dat(ii).fitDat.time) && dat(ii).fitDat.constFitChiSq < chiSq_thresh
            t1s(ii) = dat(ii).fitDat.time(1);
            a1(ii) = dat(ii).fitDat.constFitVals(2);
        end
    end
    
    d_cpd(d_cpd==0) = [];
    alph_cpd(alph_cpd==0) = [];
    
    plot(d_cpd, alph_cpd(:,1), 'ok')
    xlabel('d_{constriction} (nm)')
    ylabel('\alpha (nm^2/min)')
end

if plot_perturbed
    figure('FileName', [path '/' today '_perturbed.fig'])
    hold on
    box on
    
    for ii = 1:length(dat)
        if ~isempty(dat(ii).fitDat) && any(dat(ii).fitDat.time<dat(ii).param.t_cpd) && any(dat(ii).fitDat.time>dat(ii).param.t_cpd)
%         if ~isempty(dat(ii).fitDat) && dat(ii).fitDat.time(end)<dat(ii).param.t_cpd
%             if dat(ii).fitDat.constFitChiSq < chiSq_thresh
                plot(dat(ii).fitDat.time-dat(ii).param.t_cpd, dat(ii).fitDat.width, 'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num)], 'linew', 1)
                
                tvf = dat(ii).fitDat.time(1):0.1:dat(ii).fitDat.time(end);
                plot(tvf-dat(ii).param.t_cpd, dat(ii).fcn(dat(ii).fitDat.constFitVals,tvf), 'k', 'DisplayName', '')
%             end
        end
    end
end



