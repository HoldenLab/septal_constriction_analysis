% Author: Kevin Whitley
% Date created: 191106

% This script is used for analysing and plotting the results from multiple
% FOVs of constriction experiments.

fit_dat = 1;

plot_raw_dat = 0;
plot_fit_dat = 1;
plot_teff =  1;
plot_d_vs_alph = 0;
plot_perturbed = 0;
plot_intensity = 0;
plot_width_axial = 0;
plot_d0 = 0;

model = 'parabolic';
chiSq_thresh = 6e4;
Rsq_thresh = 0.8;

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
    ylabel('d_{septum} (nm)')
    title('Raw data')
end

if fit_dat
    
    if plot_fit_dat
        figure('FileName', [path '/' today '_fit_dat.fig'])
        hold on
        box on
    end
    
    for ii = 1:length(dat)
        
        % make sure these are both column vectors
        if size(dat(ii).cuttime,1)==1
            tfit = dat(ii).cuttime'; % old format
        else
            tfit = dat(ii).cuttime; % [min]
        end
        if size(dat(ii).cutdiams,1)==1
            dfit = dat(ii).cutdiams';
        else
            dfit = dat(ii).cutdiams; % [nm]
        end
        
        if length(tfit)<3 % too short - exclude
            continue
        end

%         if length(tfit) ~= length(dfit) % not sure why this would be the case, but happened once
%             continue
%         end
        
        % trim to remove flat ends of traces based on first derivative
        % (this is not constriction).
        w_s = smooth(dfit,10);
        t_s = smooth(tfit,10);
        
        bins = t_s(2:end)-t_s(1:end-1);
        w1d = (w_s(2:end)-w_s(1:end-1))./bins; % first derivative
        w1d_s = smooth(w1d,10);
        t1 = t_s(1:end-1) + (t_s(2:end)-t_s(1:end-1))/2;
        
        w1d_thresh = 0;
        
        halfind = floor(length(w1d_s)/2);
        indx = min([find(w1d_s(halfind:end)>=w1d_thresh, 1, 'first')+halfind length(tfit)]);
        
        % second derivative
        %         bins2 = bins(1:end-1) + (bins(2:end)-bins(1:end-1))/2;
        %         w2d = (w1d_s(2:end)-w1d_s(1:end-1))./bins2; % second derivative
        %         t2 = t1(1:end-1) + (t1(2:end)-t1(1:end-1))/2;
        
        %         w2d_thresh = 10;
        
        %         w2d_c = w2d(2:end-1);
        %         t2_c = t2(2:end-1);
        %         indx = min([find(w2d_c>=w2d_thresh, 1, 'first') length(ffit)]);
        
        tfit_c = tfit(1:indx);
        dfit_c = dfit(1:indx);
        
        % Split data into pre- and post-treatment
        
%         dat(ii).param.t_cpd = 0;% HACK!!!!!
        
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
        
        % fit pre-treatment data
        if length(dat(ii).preTimeCut) >= 5 % needs to have enough points for decent fit
            [cfit_pre, rn_pre, res, fcn] = fit_constriction_data(dat(ii).preTimeCut, dat(ii).preDiamCut, model);
            sst_pre = sum((dat(ii).preDiamCut - mean(dat(ii).preDiamCut)).^2);
            Rsq_pre = 1 - rn_pre/sst_pre;
            
            % only look at R^2 for data after constriction begins (flat
            % part at beginning really shouldn't contribute to fitting
            % error).
            onlyconpre = dat(ii).preDiamCut(dat(ii).preTimeCut >= cfit_pre(1));
            n_pre = length(onlyconpre);
            rn_onlyconpre = sum(res(dat(ii).preTimeCut >= cfit_pre(1)).^2);
            sst_pre_con = sum((onlyconpre - mean(onlyconpre)).^2);
            Rsq_pre_con = 1 - rn_onlyconpre/sst_pre_con;
        else
            cfit_pre = [];
            Rsq_pre = NaN;
            Rsq_pre_con = NaN;
            n_pre = NaN;
        end
        
        dat(ii).fitDat.preFitVals = cfit_pre;
        dat(ii).fitDat.preRsq = Rsq_pre;
        dat(ii).fitDat.preRsq_con = Rsq_pre_con;
        dat(ii).fitDat.n_pre = n_pre;
        
        % fit post-treatment data
        if length(dat(ii).postTimeCut) >= 5 && any(dat(ii).postTimeCut <= dat(ii).param.t_cpd+1) % bit shoddy, but should work for now
            [cfit_post, rn_post, res, fcn] = fit_constriction_data(dat(ii).postTimeCut, dat(ii).postDiamCut, model);
            sst_post = sum((dat(ii).postDiamCut - mean(dat(ii).postDiamCut)).^2);
            Rsq_post = 1 - rn_post/sst_post;
            
            % only look at R^2 for data after constriction begins (flat
            % part at beginning really shouldn't contribute to fitting
            % error).
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

        dat(ii).fcn = fcn;
        dat(ii).fitDat.num = ii;
        
        dat(ii).fitDat.postFitVals = cfit_post;
        dat(ii).fitDat.postRsq = Rsq_post;
        dat(ii).fitDat.postRsq_con = Rsq_post_con;
        dat(ii).fitDat.n_post = n_post;

        if plot_fit_dat

            if Rsq_pre_con > Rsq_thresh && n_pre>5
                %             hold off
                plot(dat(ii).preTime-dat(ii).param.t_cpd, dat(ii).preDiam, 'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num) '\_pre'], 'linew', 1)
                %             hold on
                tvf = dat(ii).preTime(1):0.1:dat(ii).preTime(end);
                plot(tvf-dat(ii).param.t_cpd, fcn(cfit_pre,tvf), 'k', 'DisplayName', '')
                %             teff = cfit_pre(3)^2 ./ cfit_pre(2);
                teff = d0^2 ./ cfit_pre(2);
            end
            if Rsq_post_con > Rsq_thresh && n_post>5
                plot(dat(ii).postTime-dat(ii).param.t_cpd, dat(ii).postDiam, 'DisplayName', [dat(ii).param.tracks_file(1:21) num2str(dat(ii).num) '\_post'], 'linew', 1)
                
                tvf = dat(ii).postTime(1):0.1:dat(ii).postTime(end);
                plot(tvf-dat(ii).param.t_cpd, fcn(cfit_post,tvf), 'r', 'DisplayName', '')
                %             teff = cfit(3)^2 ./ cfit(2);
                teff = d0^2 ./ cfit_post(2);
            end
        end
    end
end

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
