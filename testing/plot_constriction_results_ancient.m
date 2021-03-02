% Author: Kevin Whitley
% Date created: 191106

% This script is used for analysing and plotting the results from multiple
% FOVs of constriction experiments.

fit_dat = 1;

plot_raw_dat = 0;
plot_fit_dat = 1;
% plot_filter_dat = 0;
% plot_t_con = 0;
plot_teff =  1;
plot_d_vs_alph = 0;
plot_perturbed = 0;
plot_intensity = 0;
plot_width_axial = 0;
plot_d0 = 0;

model = 'parabolic';
chiSq_thresh = 6e4;
Rsq_thresh = 0.8;
% t0_thresh = 1000; % [min] standard error from fit
% alph_thresh = 6e11; % [nm^2/min] standard error from fit

d0 = 1000; % fix d0 for calculating teff FIGURE OUT WHAT THIS NUMBER IS!!

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

if fit_dat
    figure
    hold on
    box on
    
    rn_pre = []; rn_post = []; sst_pre = []; sst_post = [];
    for ii = 1:length(dat)
        
        ffit = dat(ii).cuttime; % [min]
        if mean(dat(ii).cutdiams)<50 % earlier data was in pixels rather than nm
            lfit = dat(ii).cutdiams * 65; % [nm]
        else
            lfit = dat(ii).cutdiams;
        end
        
        if length(dat(ii).cuttime)<3 % too short - exclude
            continue
        end

        if length(dat(ii).cuttime) ~= length(dat(ii).cutdiams) % not sure why this would be the case, but happened once
            continue
        end

%         if ~any(lfit < 400) % doesn't constrict much - fits nearly always unreliable
%             continue
%         end
        if ~any(lfit < 500) % doesn't constrict much - fits nearly always unreliable
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
        
        % second derivative
        %         bins2 = bins(1:end-1) + (bins(2:end)-bins(1:end-1))/2;
        %         w2d = (w1d_s(2:end)-w1d_s(1:end-1))./bins2; % second derivative
        %         t2 = t1(1:end-1) + (t1(2:end)-t1(1:end-1))/2;
        
        %         w2d_thresh = 10;
        
        %         w2d_c = w2d(2:end-1);
        %         t2_c = t2(2:end-1);
        %         indx = min([find(w2d_c>=w2d_thresh, 1, 'first') length(ffit)]);
        
        ffit_c = ffit(1:indx);
        lfit_c = lfit(1:indx);
        
        % Split data into pre- and post-treatment
        
%         dat(ii).param.t_cpd = 0;% HACK!!!!!
                
        dat(ii).preTime = ffit(ffit<dat(ii).param.t_cpd);
        dat(ii).preDiam = lfit(ffit<dat(ii).param.t_cpd);
        dat(ii).postTime = ffit(ffit>=dat(ii).param.t_cpd);
        dat(ii).postDiam = lfit(ffit>=dat(ii).param.t_cpd);
        
        dat(ii).preTimeCut = ffit_c(ffit_c<dat(ii).param.t_cpd);
        dat(ii).preDiamCut = lfit_c(ffit_c<dat(ii).param.t_cpd);
        dat(ii).postTimeCut = ffit_c(ffit_c>=dat(ii).param.t_cpd);
        dat(ii).postDiamCut = lfit_c(ffit_c>=dat(ii).param.t_cpd);

        if length(dat(ii).preTimeCut) >= 8 && any(dat(ii).preDiamCut < 500) % needs to have enough points for decent fit
            [cfit_pre, rn_pre(ii), ~, fcn] = fit_constriction_data(dat(ii).preTimeCut, dat(ii).preDiamCut, model);
            sst_pre(ii) = sum((dat(ii).preDiamCut - mean(dat(ii).preDiamCut)).^2);
            Rsq_pre = 1 - rn_pre(ii)/sst_pre(ii);
        else
            cfit_pre = [];
            Rsq_pre = [];
        end
%         if length(dat(ii).postTimeCut) >= 8 && any(dat(ii).postDiamCut < 500)
        if length(dat(ii).postTimeCut) >= 8 && any(dat(ii).postDiamCut < 500) && any(dat(ii).postTimeCut <= dat(ii).param.t_cpd+1) % bit shoddy, but should work for now
            [cfit_post, rn_post(ii), res, fcn] = fit_constriction_data(dat(ii).postTimeCut, dat(ii).postDiamCut, model);
            sst_post(ii) = sum((dat(ii).postDiamCut - mean(dat(ii).postDiamCut)).^2);
            Rsq_post = 1 - rn_post(ii)/sst_post(ii);
            chisq_post(ii) = rn_post(ii)/sqrt(length(dat(ii).postDiamCut));
            
            onlycon = dat(ii).postDiamCut(dat(ii).postTimeCut >= cfit_post(1));
            rn_onlycon(ii) = sum(res(dat(ii).postTimeCut >= cfit_post(1)).^2);
            sst_post_con(ii) = sum((onlycon - mean(onlycon)).^2);
            Rsq_post_con(ii) = 1 - rn_onlycon(ii)/sst_post_con(ii);
            chisq_post_con(ii) = rn_onlycon(ii)/sqrt(length(onlycon));
        else
            cfit_post = [];
            Rsq_post = [];
        end
        
        dat(ii).fcn = fcn;
        dat(ii).fitDat.num = ii;
        dat(ii).fitDat.preFitVals = cfit_pre;
        dat(ii).fitDat.postFitVals = cfit_post;
        dat(ii).fitDat.preRsq = Rsq_pre;
        dat(ii).fitDat.postRsq = Rsq_post;
        
        %         chisq = rn(ii)/sqrt(length(ffit));
        
        %         dat(ii).fitDat.constFitChiSq = chisq;
        
        if plot_fit_dat
            if Rsq_pre > Rsq_thresh
                %             hold off
                plot(dat(ii).preTime-dat(ii).param.t_cpd, dat(ii).preDiam, 'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num) '\_pre'], 'linew', 1)
                %             hold on
                tvf = dat(ii).preTime(1):0.1:dat(ii).preTime(end);
                plot(tvf-dat(ii).param.t_cpd, fcn(cfit_pre,tvf), 'k', 'DisplayName', '')
                %             teff = cfit_pre(3)^2 ./ cfit_pre(2);
                teff = d0^2 ./ cfit_pre(2);
            end
            if Rsq_post > Rsq_thresh
                plot(dat(ii).postTime-dat(ii).param.t_cpd, dat(ii).postDiam, 'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num) '\_post'], 'linew', 1)
                
                tvf = dat(ii).postTime(1):0.1:dat(ii).postTime(end);
                plot(tvf-dat(ii).param.t_cpd, fcn(cfit_post,tvf), 'r', 'DisplayName', '')
                %             teff = cfit(3)^2 ./ cfit(2);
                teff = d0^2 ./ cfit_post(2);
            end
        end
    end
end

% Deprecated 200216
% if plot_fit_dat
%     figure
%     hold on
%     box on
%     
%     %     fov = alldat.fov(jj);
%     for ii = 1:length(dat)
%         plot(dat.fitDat(ii).time, dat.fitDat(ii).width, 'DisplayName', num2str(dat.fitDat(ii).num))
%         
%         tvf = fov.fitDat(ii).time(1):0.1:fov.fitDat(ii).time(end);
%         plot(tvf, fov.fcn(fov.fitDat(ii).constFitVals,tvf), 'k', 'DisplayName', '')
%     end
%     tv = 0:1200;
%     plot(ones(length(tv),1)*fov.param.t_cpd, tv, 'r')
%     
%     xlabel('Time (min)')
%     ylabel('d_{septum} (nm)')
%     title(['Fitted, unfiltered data for ' strrep(fov.fname,'_','\_')])
% end
% 
% if plot_filter_dat
%     figure
%     hold on
%     box on
%     
%     for jj = 1:length(alldat.fov)
%         %         fov = alldat.fov(jj);
%         fov = alldat2(jj);
%         for ii = 1:length(fov.fitDat)
%             if fov.fitDat(ii).constFitChiSq < chiSq_thresh
%                 plot(fov.fitDat(ii).time, fov.fitDat(ii).width, 'DisplayName', num2str(fov.fitDat(ii).num), 'linew', 1)
%                 
%                 tvf = fov.fitDat(ii).time(1):0.1:fov.fitDat(ii).time(end);
%                 plot(tvf, fov.fcn(fov.fitDat(ii).constFitVals,tvf), 'k', 'DisplayName', '')
%             end
%         end
%     end
%     tv = 0:1200;
%     plot(ones(length(tv),1)*fov.param.t_cpd, tv, 'r')
%     
%     xlabel('Time (min)')
%     ylabel('d_{septum} (nm)')
%     title(['Fitted, filtered data for ' strrep(fov(1).param.tracks_file,'_','\_')])
% end
% 
% if plot_t_con
%     figure('FileName', [path '/' today '_t_con.fig'])
%     hold on
%     box on
%     
%     t_con = [];
%     for jj = 1:length(alldat.fov)
%         fov = alldat.fov(jj);
%         for ii = 1:length(fov.fitDat)
%             if fov.fitDat(ii).constFitChiSq < chiSq_thresh
%                 t_con = [t_con; fov.fitDat(ii).time(end) - fov.fitDat(ii).constFitVals(1) fov.fitDat(ii).time(end)];
%             end
%         end
%     end
%     
%     t_con(t_con(:,1)<1,:) = [];
%     
%     t_con_pre = t_con(t_con(:,2)<fov.param.t_cpd,1);
%     t_con_post = t_con(t_con(:,2)>fov.param.t_cpd,1);
%     tpr = histogram(t_con_pre,5:2:40);
%     tpo = histogram(t_con_post,5:2:40);
%     set(tpr,'facea',0.5)
%     set(tpo,'facea',0.5)
%     xlabel('t_{constriction} (min)')
%     ylabel('Counts')
% end

if plot_teff
    figure('FileName', [path '/' today '_teff.fig'])
    hold on
    box on
    
    teff_pre = []; teff_post = [];
    for ii = 1:length(dat)
%         if ~isempty(dat(ii).fitDat) && dat(ii).fitDat.Rsq > Rsq_thresh
%             teff2 = dat(ii).fitDat.constFitVals(3)^2 ./ dat(ii).fitDat.constFitVals(2);
%             teff = [teff; dat(ii).fitDat.constFitVals(2) dat(ii).fitDat.time(end) dat(ii).param.t_cpd teff2];
%         end
        if ~isempty(dat(ii).fitDat) && ~isempty(dat(ii).fitDat.preRsq) && dat(ii).fitDat.preRsq > Rsq_thresh
%             teff_pre(ii) = dat(ii).fitDat.preFitVals(3)^2 ./ dat(ii).fitDat.preFitVals(2);
            teff_pre(ii) = d0^2 ./ dat(ii).fitDat.preFitVals(2);
        else
            teff_pre(ii) = NaN;
        end
        if ~isempty(dat(ii).fitDat) && ~isempty(dat(ii).fitDat.postRsq) && dat(ii).fitDat.postRsq > Rsq_thresh
%             teff_post(ii) = dat(ii).fitDat.postFitVals(3)^2 ./ dat(ii).fitDat.postFitVals(2);
            teff_post(ii) = d0^2 ./ dat(ii).fitDat.postFitVals(2);
        else
            teff_post(ii) = NaN;
        end
    end
    
%     teff_pre = teff(teff(:,2)<teff(:,3),4);
%     teff_post = teff(teff(:,2)>teff(:,3),4);
    %     alph_pre = teff(teff(:,2)<teff(:,3),1);
    %     alph_post = teff(teff(:,2)>teff(:,3),1);
    
    %     apr = histogram(alph_pre,1e4:5e3:1e5);
    %     apo = histogram(alph_post,1e4:5e3:1e5);
    %     set(apr,'facea',0.5)
    %     set(apo,'facea',0.5)
    %     xlabel('\alpha (nm^2/min)')
    %     ylabel('Counts')
    
    %     alph_pre = [alph_pre; ones(length(alph_post)-length(alph_pre),1)*NaN];
    %     alph_post = [alph_post; ones(length(alph_pre)-length(alph_post),1)*NaN];
%     teff_pre = [teff_pre; ones(length(teff_post)-length(teff_pre),1)*NaN];
%     teff_post = [teff_post; ones(length(teff_pre)-length(teff_post),1)*NaN];
    
    if any(~isnan(teff_pre))
        violinplot([teff_pre' teff_post'])
    else
        violinplot(teff_post)
    end
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

if plot_intensity
    figure('FileName', [path '/' today '_intensity.fig'])
    hold on
    box on
    
    for ii = 1:length(dat)
        if any(dat(ii).cuttime<dat(ii).param.t_cpd) && any(dat(ii).cuttime>dat(ii).param.t_cpd)
            plot(dat(ii).cuttime, dat(ii).cutint, ...
                'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num)])
        end
    end
end

if plot_width_axial
    figure('FileName', [path '/' today '_width_axial.fig'])
    hold on
    box on
    
    for ii = 1:length(dat)
        plot(dat(ii).cuttime_ax, dat(ii).cutdiams_ax, ...
            'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num)])
    end
end

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
