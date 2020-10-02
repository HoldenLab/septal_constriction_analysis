% Author: Kevin Whitley
% Date created: 191203

% This script takes a dataset variable 'results2' containing information on
% septal width, thickness, intensity, and whether the septa
% condensed/constricted or not post-treatment, and plots various things
% from it (see options).

% The matrix results2 should be in format: 1-date, 2-file number, 3-FOV,
% 4-ring number, 5-center of line profile [pix], 6-radius [pix], 7-fit
% error, 8-constricted?, 9-nascent ring radius [pix], 10-axial width FWHM
% [pix], 11-axial fit error, 12-median pixel intensity of HT-PBP2B,
% 13-manual classification of state, 14-median pixel intensity of FtsZ-GFP,
% 15-power of laser illuminating PBP2B (uW), 16-power of laser illuminating
% FtsZ (uW) 17-condensed? 18-low 2B signal increased?

dat = results2;
ud = [];

% Analysis/plotting options

% make bar plots of % septa constricted post-treatment vs. various things
plot_f_con_vs_diameter = 0; % vs. septal diameter
plot_f_con_vs_thickness = 0; % vs. septal thickness
plot_f_con_vs_2Bintensity = 0; % vs. PBP2B intensity
plot_f_con_vs_Zintensity = 0; % vs. FtsZ intensity
plot_f_con_manualclass = 0; % vs. division stage, manually classified (column 13)
plot_f_con_autoclass = 0; % vs. division stage, automatically classified by diameter and thickness
plot_f_con_2Bint_norm = 1; % vs. low/high PBP2B intensity (normalized if possible)
plot_f_con_Zint_norm = 0; % vs. low/high FtsZ intensity (normalized if possible)

% make scatter plots of various properties and whether septa continued
% division post-treatment
plot_scatter_diameter_thickness = 0; % septal diameter vs. thickness
plot_scatter_diameter_2Bintensity = 0; % septal diameter vs. PBP2B intensity (normalized if possible)
plot_scatter_diameter_Zintensity = 0; % septal diameter vs. FtsZ intensity (normalized if possible)
plot_scatter_2Bintensity_Zintensity = 0; % PBP2B intensity vs. FtsZ intensity (both normalized if possible)

plot_f_condense_thickness = 0; % make stacked bar plot of % septa condensed post-treatment vs. septal thickness
plot_scatter_diameter_thickness_condense = 0; % make scatter plot of septal diameter vs. thickness and whether septa condensed post-treatment

plot_f_low2B_increase = 0; % make stacked bar plot of % septa where low PBP2B signal increased post-treatment

ud.use_relative_diameter = 1; % use relative septal diameter (otherwise use absolute)
ud.is_2color = 1; % is this the 2-color strain?

if ud.is_2color
    thickness_thresh = 450; % [nm] 2-color (poor SNR)
else
    thickness_thresh = 400; % [nm] 1-color
end

% Thresholds for data filters
cent = median(dat(~isnan(dat(:,5)),5));
sig_cent = nanstd(dat(:,5));
ud.cent_cut = [cent-2*sig_cent cent+2*sig_cent]; % [pix] cut on fitted center of line profile (if off-center, likely did not fit whole ring)
% ud.diam_cut = [700 1200] /65/2; % [pix] cut on radial diameter
ud.diam_cut = [0 Inf] /65/2; % [pix] cut on radial diameter
ud.diam_err_cut = [0 0.7]; % cut on radial fitting error
% ud.thick_cut = [0 thickness_thresh] /65; % [pix] cut on axial width
ud.thick_cut = [0 Inf] /65; % [pix] cut on axial width
ud.thick_err_cut = [0 1]; % cut on error in fitting for axial diameters
ud.int_2b_cut = [0 Inf]; % cut on 2B intensity
ud.class_cut = [0 Inf]; % cut on classification
ud.int_Z_cut = [0 Inf]; % cut on Z intensity
if ud.use_relative_diameter
    radcut = 0.9;
    ud.rel_diam_cut = [0 1.3]; % cut on relative diameter
else
    radcut = 700;
end

% Add relative septal diameters to data matrix based on each individual d0
% (column 9, if available)
if ud.use_relative_diameter
    lastcol = size(dat,2) + 1; % add a new column to data matrix to put in relative diameters
    dat(:,lastcol) = dat(:,6) ./ dat(:,9); % relative diameters
    dat(dat(:,lastcol)==Inf,lastcol) = 1; % If nascent ring diameter=0, must be a nascent ring. this is max diameter.
end

% replace intensities with normalized intensities (per uW excitation power)
if size(dat,2)>15 && any(dat(:,16))
    dat(:,14) = dat(:,14) ./ dat(:,16); % normalized intensities of FtsZ
end
if size(dat,2)>14 && any(dat(:,15))
    dat(:,12) = dat(:,12) ./ dat(:,15); % normalized intensities of PBP2B
end

% Filter data
dat = cutmat(dat, 5, ud.cent_cut); % cut on fitted center position
dat = cutmat(dat, 6, ud.diam_cut); % cut on diameter
dat = cutmat(dat, 7, ud.diam_err_cut); % cut on error in diameter fit
dat(isnan(dat(:,8)),:) = []; % remove rows where can't tell if septum constricted (e.g. moved out of FOV)
dat = cutmat(dat, 10, ud.thick_cut); % cut on thickness
dat = cutmat(dat, 11, ud.thick_err_cut); % cut on error in thickness fit
dat = cutmat(dat, 12, ud.int_2b_cut); % cut on PBP2B intensity
dat = cutmat(dat, 13, ud.class_cut); % cut on manual classification
dat = cutmat(dat, 14, ud.int_Z_cut); % cut on FtsZ intensity
if ud.use_relative_diameter
    if size(dat,2)>18
        dat = cutmat(dat, 19, ud.cent_cut); % cut on nascent ring center position
        dat = cutmat(dat, 20, ud.diam_err_cut); % cut on nascent ring diameter fit error
    end
    dat = cutmat(dat, lastcol, ud.rel_diam_cut); % cut on relative diameter
end

if plot_f_con_vs_diameter
    
    if ud.use_relative_diameter
        drange = [0.01 1.3]; % relative
        dbsize = 0.1;
    else
        drange = [0 1500]; % absolute
        dbsize = 100;
    end
    dbins = drange(1):dbsize:drange(end); % diameter bins
    
    fcon = []; se_fcon = []; n = [];
    for ii = 1:length(dbins)

        if ud.use_relative_diameter
            datind = dbins(ii)-dbsize/2 < dat(:,lastcol) & dat(:,lastcol) < dbins(ii)+dbsize/2; % relative
        else
            datind = dbins(ii)-dbsize/2 < dat(:,6)*2*65 & dat(:,6)*2*65 < dbins(ii)+dbsize/2; % absolute
        end
        
        cutdat = dat(datind,:);
        
        if length(cutdat(:,8)) > 1
            b_fcon = bootstrp(1000, @(x) sum(x)./length(x), cutdat(:,8));
            fcon(ii) = mean(b_fcon);
            se_fcon(ii) = std(b_fcon);
            n(ii) = length(cutdat(:,8));
        elseif ~isempty(cutdat)
            fcon(ii) = cutdat(:,8);
            se_fcon(ii) = NaN;
        else
            fcon(ii) = NaN;
            se_fcon(ii) = NaN;
        end
        
    end
    
    figure
    hold on
    box on
    
    bar(dbins, fcon*100)
    errorbar(dbins, fcon*100, se_fcon*100/2, 'k', 'LineStyle', 'none')
%     set(gca,'XDir','reverse')
    ylim([0 105])
    if ud.use_relative_diameter
        xlabel('Relative septal diameter')
    else
        xlabel('Septal diameter (nm)')
    end
    ylabel('% Continuing constriction')
    set(gcf, 'UserData', ud)
end

if plot_f_con_vs_thickness
    
    arange = [0 1200];
    absize = 100;
    
    abins = arange(1):absize:arange(end); % [nm] diameter bins
    
    fcon = []; se_fcon = []; n = [];
    for ii = 1:length(abins)
        
        datind = abins(ii)-absize/2 < dat(:,10)*65 & dat(:,10)*65 < abins(ii)+absize/2; % relative
        cutdat = dat(datind,:);
        
        if length(cutdat(:,8)) > 1
            b_fcon = bootstrp(1000, @(x) sum(x)./length(x), cutdat(:,8));
            fcon(ii) = mean(b_fcon);
            se_fcon(ii) = std(b_fcon);
            n(ii) = length(cutdat(:,8));
        elseif ~isempty(cutdat)
            fcon(ii) = cutdat(:,8);
            se_fcon(ii) = NaN;
        else
            fcon(ii) = NaN;
            se_fcon(ii) = NaN;
        end
        
    end
    
    figure
    hold on
    box on
    
    bar(abins, fcon*100)
    errorbar(abins, fcon*100, se_fcon*100/2, 'k', 'LineStyle','none')
    set(gca,'XDir','reverse')
    ylim([0 105])
    xlabel('Septal thickness (nm)')
    ylabel('% Continuing constriction')
    set(gcf, 'UserData', ud)
end

if plot_f_con_vs_2Bintensity
    
    irange = [min(dat(:,12)) max(dat(:,12))];
    ibsize = (irange(2)-irange(1))/10;
    
    ibins = irange(1):ibsize:irange(end); % [nm] diameter bins
    
    fcon = []; se_fcon = []; n = [];
    for ii = 1:length(ibins)
        
        datind = ibins(ii)-ibsize/2 < dat(:,12) & dat(:,12) < ibins(ii)+ibsize/2; % relative
        cutdat = dat(datind,:);
        
        if length(cutdat(:,8)) > 1
            b_fcon = bootstrp(1000, @(x) sum(x)./length(x), cutdat(:,8));
            fcon(ii) = mean(b_fcon);
            se_fcon(ii) = std(b_fcon);
            n(ii) = length(cutdat(:,8));
        elseif ~isempty(cutdat)
            fcon(ii) = cutdat(:,8);
            se_fcon(ii) = NaN;
        else
            fcon(ii) = NaN;
            se_fcon(ii) = NaN;
        end
        
    end
    
    figure
    hold on
    box on
    
    bar(ibins, fcon*100)
    errorbar(ibins, fcon*100, se_fcon*100/2, 'k', 'LineStyle', 'none')
    set(gca,'XDir','reverse')
    ylim([0 105])
    xlabel('Median PBP2B pixel intensity')
    ylabel('% Continuing constriction')
    title('% of cells continuing constriction post-PC19 vs. PBP2B intensity.')
    set(gcf, 'UserData', ud)
end

if plot_f_con_vs_Zintensity
    
    irange = [0 max(dat(:,14))];
    ibsize = (irange(2)-irange(1))/10;
    
    ibins = irange(1):ibsize:irange(end); % [nm] diameter bins
    
    fcon = []; se_fcon = []; n = [];
    for ii = 1:length(ibins)
        
        datind = ibins(ii)-ibsize/2 < dat(:,14) & dat(:,14) < ibins(ii)+ibsize/2; % relative
        cutdat = dat(datind,:);
        
        if length(cutdat(:,8)) > 1
            b_fcon = bootstrp(1000, @(x) sum(x)./length(x), cutdat(:,8));
            fcon(ii) = mean(b_fcon);
            se_fcon(ii) = std(b_fcon);
            n(ii) = length(cutdat(:,8));
        else
            fcon(ii) = NaN;
            se_fcon(ii) = NaN;
        end
        
    end
    
    figure
    hold on
    box on
    
    bar(ibins, fcon*100)
    errorbar(ibins, fcon*100, se_fcon*100/2, 'k', 'LineStyle', 'none')
    set(gca,'XDir','reverse')
    ylim([0 105])
    xlabel('Median FtsZ-GFP pixel intensity')
    ylabel('% Continuing constriction')
    title('% of cells continuing constriction post-PC19 vs. FtsZ-GFP intensity.')
    set(gcf, 'UserData', ud)
end

if plot_f_con_manualclass
    
    mans = dat(:,13);
    
    fcon = []; se_fcon = []; n = [];
    for ii = 0:2
        
        cutdat = dat(mans==ii,:);
        
        if length(cutdat(:,8)) > 1
            b_fcon = bootstrp(1000, @(x) sum(x)./length(x), cutdat(:,8));
            fcon(ii+1) = mean(b_fcon);
            se_fcon(ii+1) = std(b_fcon);
            n(ii+1) = length(cutdat(:,8));
        else
            fcon(ii+1) = NaN;
            se_fcon(ii+1) = NaN;
        end
    end
    
    figure
    hold on
    box on
    
    bar(0:2, fcon*100)
    errorbar(0:2, fcon*100, se_fcon*100/2, 'k', 'LineStyle','none')
    ylim([0 105])
    ylabel('% Continuing constriction')
    set(gcf, 'UserData', ud)
end

if plot_f_con_autoclass
    
    axcut_nasc = [thickness_thresh 1e5] /65;
    nasc = cutmat(dat, 10, axcut_nasc); % matrix of nascent rings
    
    axcut_mat = [0 thickness_thresh] /65;
    matcon = cutmat(dat, 10, axcut_mat); % matrix of mature and constricting rings
    
    if ud.use_relative_diameter
        radcut_mat = [radcut 1.3];        
        radcut_cons = [0 radcut];
        cutcol = lastcol;
    else
        radcut_mat = [radcut 1200] /65/2;
        radcut_cons = [0 radcut] /65/2;
        cutcol = 6;
    end
    
    mats = cutmat(matcon, cutcol, radcut_mat); % matrix of mature rings
    cons = cutmat(matcon, cutcol, radcut_cons); % matrix of constricting rings
    
    fcon=[]; se_fcon=[];
    if length(nasc(:,8)) > 1
%         b_nasc = bootstrp(1000, @(x) sum(x)./length(x), nasc(:,8));
%         fcon(1) = mean(b_nasc);
%         se_fcon(1) = std(b_nasc);
        n(1) = length(nasc(:,8));
        
        fcon_nasc = sum(nasc(:,8))/n(1);
        fcon(1,:) = [fcon_nasc 1-fcon_nasc];
    elseif ~isempty(nasc)
        fcon(1,:) = [1 0];
%         se_fcon(1) = NaN;
    else
        fcon(1,:) = [NaN NaN];
%         se_fcon(1) = NaN;
    end

    if length(mats(:,8)) > 1
%         b_mats = bootstrp(1000, @(x) sum(x)./length(x), mats(:,8));
%         fcon(2) = mean(b_mats);
%         se_fcon(2) = std(b_mats);
        n(2) = length(mats(:,8));
        
        fcon_mat = sum(mats(:,8))/n(2);
        fcon(2,:) = [fcon_mat 1-fcon_mat];
    elseif ~isempty(mats)
        fcon(2,:) = [1 0];
%         se_fcon(2) = NaN;
    else
        fcon(2,:) = NaN;
%         se_fcon(2) = NaN;
    end

    if length(cons(:,8)) > 1
%         b_cons = bootstrp(1000, @(x) sum(x)./length(x), cons(:,8));
%         fcon(3) = mean(b_cons);
%         se_fcon(3) = std(b_cons);
        n(3) = length(cons(:,8));
        
        fcon_con = sum(cons(:,8))/n(3);
        fcon(3,:) = [fcon_con 1-fcon_con];
    elseif ~isempty(cons)
        fcon(3,:) = [1 0];
%         se_fcon(3) = NaN;
    else
        fcon(3,:) = [NaN NaN];
%         se_fcon(3) = NaN;
    end
    
    figure
    hold on
    box on
    
    bar(fcon*100, 'stacked')
%     bar(0:2, fcon*100)
%     errorbar(0:2, fcon*100, se_fcon*100/2, 'k', 'LineStyle','none')
    ylim([0 105])
    ylabel('% Continuing constriction')
    set(gcf, 'UserData', ud)
end

if plot_f_con_2Bint_norm

    if ud.is_2color
        int_cut = 0.15; % 2-color (poor SNR)
    else
        int_cut = 0.1; % 1-color
    end
    
    low_range = [0 int_cut]; % low 2B signal
    cutdat_low = cutmat(dat, 12, low_range);
    
    high_range = [int_cut 1e5]; % high 2B signal
    cutdat_high = cutmat(dat, 12, high_range);
    
    fcon=[]; se_fcon=[]; n=[];
    if length(cutdat_low(:,8)) > 1
%         b_no = bootstrp(1000, @(x) sum(x)./length(x), cutdat_no(:,8));
%         fcon(1) = mean(b_no);
%         se_fcon(1) = std(b_no);
        n(1) = length(cutdat_low(:,8));
        
        fcon_low = sum(cutdat_low(:,8))/n(1);
        fcon(1,:) = [fcon_low 1-fcon_low];
    elseif ~isempty(cutdat_low)
        fcon(1,:) = [0 1];
%         fcon(1) = cutdat_no(:,8);
        se_fcon(1) = NaN;
    else
        fcon(1,:) = [NaN NaN];
        se_fcon(1) = NaN;
    end
    
    if length(cutdat_high(:,8)) > 1
%         b_yes = bootstrp(1000, @(x) sum(x)./length(x), cutdat_yes(:,8));
%         fcon(2) = mean(b_yes);
%         se_fcon(2) = std(b_yes);
        n(2) = length(cutdat_high(:,8));
        
        fcon_high = sum(cutdat_high(:,8))/n(2);
        fcon(2,:) = [fcon_high 1-fcon_high];
    elseif ~isempty(cutdat_high)
        fcon(2,:) = [0 1];
%         se_fcon(2) = NaN;
    else
        fcon(2,:) = [NaN NaN];
%         se_fcon(2) = NaN;
    end
    
    figure
    hold on
    box on
    
%     bar(0:1, fcon*100)
    bar(fcon*100, 'stacked')
%     errorbar(0:1, fcon*100, se_fcon*100/2, 'k', 'LineStyle','none')
    ylim([0 105])
    ylabel('% Continuing constriction')
    set(gcf, 'UserData', ud)
end

if plot_f_con_Zint_norm

    low_range = [0 1.2]; % low Z signal
    cutdat_low = cutmat(dat, 14, low_range);
    
    high_range = [1.2 1e5]; % high Z signal
    cutdat_high = cutmat(dat, 14, high_range);
    
    fcon=[]; se_fcon=[];
    if length(cutdat_low(:,8)) > 1
%         b_no = bootstrp(1000, @(x) sum(x)./length(x), cutdat_no(:,8));
%         fcon(1) = mean(b_no);
%         se_fcon(1) = std(b_no);
        n(1) = length(cutdat_low(:,8));
        
        fcon_low = sum(cutdat_low(:,8))/n(1);
        fcon(1,:) = [fcon_low 1-fcon_low];
    elseif ~isempty(cutdat_low)
        fcon(1,:) = [0 1];
%         fcon(1) = cutdat_no(:,8);
        se_fcon(1) = NaN;
    else
        fcon(1,:) = [NaN NaN];
        se_fcon(1) = NaN;
    end
    
    if length(cutdat_high(:,8)) > 1
%         b_yes = bootstrp(1000, @(x) sum(x)./length(x), cutdat_yes(:,8));
%         fcon(2) = mean(b_yes);
%         se_fcon(2) = std(b_yes);
        n(2) = length(cutdat_high(:,8));
        
        fcon_high = sum(cutdat_high(:,8))/n(2);
        fcon(2,:) = [fcon_high 1-fcon_high];
    elseif ~isempty(cutdat_high)
        fcon(2,:) = [0 1];
%         se_fcon(2) = NaN;
    else
        fcon(2,:) = [NaN NaN];
%         se_fcon(2) = NaN;
    end
    
    figure
    hold on
    box on
    
%     bar(0:1, fcon*100)
    bar(fcon*100, 'stacked')
%     errorbar(0:1, fcon*100, se_fcon*100/2, 'k', 'LineStyle','none')
    ylim([0 105])
    ylabel('% Continuing constriction')
    set(gcf, 'UserData', ud)
end

if plot_scatter_diameter_thickness
    
    yups = dat(dat(:,8)==1,:); % continued division
    nopes = dat(dat(:,8)==0,:); % did not continue division
    
    figure
    hold on
    box on
    
    if ud.use_relative_diameter
        scatter(yups(:,lastcol), yups(:,10)*65, 'MarkerEdgeColor', 'k',...
            'MarkerFaceColor', 'b', 'MarkerFaceAlpha', 0.5, 'SizeData', 50)
        scatter(nopes(:,lastcol), nopes(:,10)*65, 'MarkerEdgeColor', 'k',...
            'MarkerFaceColor', 'r', 'MarkerFaceAlpha', 0.5, 'SizeData', 50)
    else
        scatter(yups(:,6)*65*2, yups(:,10)*65, 'MarkerEdgeColor', 'k',...
            'MarkerFaceColor', 'b', 'MarkerFaceAlpha', 0.5, 'SizeData', 50)
        scatter(nopes(:,6)*65*2, nopes(:,10)*65, 'MarkerEdgeColor', 'k',...
            'MarkerFaceColor', 'r', 'MarkerFaceAlpha', 0.5, 'SizeData', 50)
    end

    legend('Continued constricting', 'Did not constrict')
    if ud.use_relative_diameter
        xlabel('Relative septal diameter')
    else
        xlabel('Diameter (nm)')
    end
    ylabel('Ring thickness (nm)')
    title('Thickness and diameter measurements pre-PC19 along with post-PC19 determination of constriction fate.')
    set(gcf, 'UserData', ud)
end

if plot_scatter_diameter_2Bintensity

    yups = dat(dat(:,8)==1,:); % continued division
    nopes = dat(dat(:,8)==0,:); % did not continue division
    
    figure
    hold on
    box on
    
    if ud.use_relative_diameter
        scatter(yups(:,lastcol), yups(:,12), 'MarkerEdgeColor', 'k',...
            'MarkerFaceColor', 'b', 'MarkerFaceAlpha', 0.5, 'SizeData', 50)
        scatter(nopes(:,lastcol), nopes(:,12), 'MarkerEdgeColor', 'k',...
            'MarkerFaceColor', 'r', 'MarkerFaceAlpha', 0.5, 'SizeData', 50)
    else
        scatter(yups(:,6)*65*2, yups(:,12), 'MarkerEdgeColor', 'k',...
            'MarkerFaceColor', 'b', 'MarkerFaceAlpha', 0.5, 'SizeData', 50)
        scatter(nopes(:,6)*65*2, nopes(:,12), 'MarkerEdgeColor', 'k',...
            'MarkerFaceColor', 'r', 'MarkerFaceAlpha', 0.5, 'SizeData', 50)
    end
    
    legend('Continued constricting', 'Did not constrict')
    xlabel('Diameter (nm)')
    ylabel('PBP2B median pixel intensity')
    title('Diameter and 2B intensity measurements pre-PC19 along with post-PC19 determination of division fate.')
    set(gcf, 'UserData', ud)
end

if plot_scatter_diameter_Zintensity
    
    yups = dat(dat(:,8)==1,:); % continued division
    nopes = dat(dat(:,8)==0,:); % did not continue division
    
    figure
    hold on
    box on
    
    if ud.use_relative_diameter
        scatter(yups(:,lastcol), yups(:,14), 'MarkerEdgeColor', 'k',...
            'MarkerFaceColor', 'b', 'MarkerFaceAlpha', 0.5, 'SizeData', 50)
        scatter(nopes(:,lastcol), nopes(:,14), 'MarkerEdgeColor', 'k',...
            'MarkerFaceColor', 'r', 'MarkerFaceAlpha', 0.5, 'SizeData', 50)
    else
        scatter(yups(:,6)*65*2, yups(:,14), 'MarkerEdgeColor', 'k',...
            'MarkerFaceColor', 'b', 'MarkerFaceAlpha', 0.5, 'SizeData', 50)
        scatter(nopes(:,6)*65*2, nopes(:,14), 'MarkerEdgeColor', 'k',...
            'MarkerFaceColor', 'r', 'MarkerFaceAlpha', 0.5, 'SizeData', 50)
    end
    
    legend('Continued constricting', 'Did not constrict')
    xlabel('Diameter (nm)')
    ylabel('FtsZ-GFP median pixel intensity')
    title('Diameter and Z intensity measurements pre-PC19 along with post-PC19 determination of division fate.')
    set(gcf, 'UserData', ud)
end

if plot_scatter_2Bintensity_Zintensity
    
    yups = dat(dat(:,8)==1,:); % continued division
    nopes = dat(dat(:,8)==0,:); % did not continue division
    
    figure
    hold on
    box on
    
    scatter(yups(:,12), yups(:,14), 'MarkerEdgeColor', 'k',...
        'MarkerFaceColor', 'b', 'MarkerFaceAlpha', 0.5, 'SizeData', 50)
    scatter(nopes(:,12), nopes(:,14), 'MarkerEdgeColor', 'k',...
        'MarkerFaceColor', 'r', 'MarkerFaceAlpha', 0.5, 'SizeData', 50)
    
    legend('Continued constricting', 'Did not constrict')
    xlabel('PBP2B median pixel intensity')
    ylabel('FtsZ-GFP median pixel intensity')
    title('PBP2B and FtsZ-GFP intensity measurements pre-PC19 along with post-PC19 determination of division fate.')
    set(gcf, 'UserData', ud)
end

if plot_f_condense_thickness
    
    thick_range = [thickness_thresh 1e5]/65; % thick
    cutdat_thick = cutmat(dat, 10, thick_range);
    cutdat_thick(isnan(cutdat_thick(:,17)),:) = [];
    
    thin_range = [0 thickness_thresh]/65; % thin
    cutdat_thin = cutmat(dat, 10, thin_range);
    cutdat_thin(isnan(cutdat_thin(:,17)),:) = [];
    
    fcon=[]; se_fcon=[];
    if length(cutdat_thick(:,17)) > 1
%         b_thick = bootstrp(1000, @(x) sum(x)./length(x), cutdat_thick(:,17));
%         fcon(1) = nanmean(b_thick);
%         se_fcon(1) = nanstd(b_thick);
        n(1) = length(cutdat_thick(~isnan(cutdat_thick(:,17)),17));
        
        fcon_thick = sum(cutdat_thick(:,17))/n(1);
        fcon(1,:) = [fcon_thick 1-fcon_thick];
    elseif ~isempty(cutdat_thick)
        fcon(1,:) = [0 1];
%         se_fcon(1) = NaN;
    else
        fcon(1,:) = [NaN NaN];
%         se_fcon(1) = NaN;
    end
    
    if length(cutdat_thin(:,17)) > 1
%         b_thin = bootstrp(1000, @(x) sum(x)./length(x), cutdat_thin(:,17));
%         fcon(2) = nanmean(b_thin);
%         se_fcon(2) = nanstd(b_thin);
        n(2) = length(cutdat_thin(~isnan(cutdat_thin(:,17)),17));
        
        fcon_thin = sum(cutdat_thin(:,17))/n(2);
        fcon(2,:) = [fcon_thin 1-fcon_thin];
    elseif ~isempty(cutdat_thin)
        fcon(2,:) = [0 1];
%         se_fcon(2) = NaN;
    else
        fcon(2) = [NaN NaN];
%         se_fcon(2) = NaN;
    end
    
    figure
    hold on
    box on
    
%     bar(0:1, fcon*100)
    bar(fcon*100, 'stacked')
%     errorbar(0:1, fcon*100, se_fcon*100/2, 'k', 'LineStyle','none')
    ylim([0 105])
    ylabel('% Condensed')
    set(gcf, 'UserData', ud)
end

if plot_scatter_diameter_thickness_condense
    
    yups = dat(dat(:,17)==1,:); % condensed
    nopes = dat(dat(:,17)==0,:); % did not condense
    
    figure
    hold on
    box on
    
    if ud.use_relative_diameter
        scatter(yups(:,17), yups(:,10)*65, 'MarkerEdgeColor', 'k',...
            'MarkerFaceColor', 'b', 'MarkerFaceAlpha', 0.5, 'SizeData', 50)
        scatter(nopes(:,17), nopes(:,10)*65, 'MarkerEdgeColor', 'k',...
            'MarkerFaceColor', 'r', 'MarkerFaceAlpha', 0.5, 'SizeData', 50)
    else
        scatter(yups(:,6)*65*2, yups(:,10)*65, 'MarkerEdgeColor', 'k',...
            'MarkerFaceColor', 'b', 'MarkerFaceAlpha', 0.5, 'SizeData', 50)
        scatter(nopes(:,6)*65*2, nopes(:,10)*65, 'MarkerEdgeColor', 'k',...
            'MarkerFaceColor', 'r', 'MarkerFaceAlpha', 0.5, 'SizeData', 50)
    end

    legend('Condensed', 'Did not condense')
    if ud.use_relative_diameter
        xlabel('Relative septal diameter')
    else
        xlabel('Diameter (nm)')
    end
    ylabel('Ring thickness (nm)')
    title('Thickness and diameter measurements pre-PC19 along with post-PC19 determination of condensation fate.')
    set(gcf, 'UserData', ud)
end

if plot_f_low2B_increase

    thick_range = [thickness_thresh 1e5]/65; % thick
    cutdat_thick = cutmat(dat, 10, thick_range);
    cutdat_thick(isnan(cutdat_thick(:,18)),:) = [];
    
    thin_range = [0 thickness_thresh]/65; % thin
    cutdat_thin = cutmat(dat, 10, thin_range);
    cutdat_thin(isnan(cutdat_thin(:,18)),:) = [];
    
    fcon=[]; se_fcon=[];
    if length(cutdat_thick(:,18)) > 1
%         b_no = bootstrp(1000, @(x) sum(x)./length(x), cutdat_thick(:,17));
%         fcon(1) = nanmean(b_thick);
%         se_fcon(1) = nanstd(b_thick);
        n(1) = length(cutdat_thick(~isnan(cutdat_thick(:,18)),18));
        
        fcon_thick = sum(cutdat_thick(:,18))/n(1);
        fcon(1,:) = [fcon_thick 1-fcon_thick];
    elseif ~isempty(cutdat_thick)
        fcon(1,:) = [0 1];
%         se_fcon(1) = NaN;
    else
        fcon(1,:) = [NaN NaN];
%         se_fcon(1) = NaN;
    end
    
    if length(cutdat_thin(:,18)) > 1
%         b_thin = bootstrp(1000, @(x) sum(x)./length(x), cutdat_thin(:,17));
%         fcon(2) = nanmean(b_thin);
%         se_fcon(2) = nanstd(b_thin);
        n(2) = length(cutdat_thin(~isnan(cutdat_thin(:,18)),18));
        
        fcon_thin = sum(cutdat_thin(:,18))/n(2);
        fcon(2,:) = [fcon_thin 1-fcon_thin];
    elseif ~isempty(cutdat_thin)
        fcon(2,:) = [0 1];
%         se_fcon(2) = NaN;
    else
        fcon(2) = [NaN NaN];
%         se_fcon(2) = NaN;
    end
    
    figure
    hold on
    box on
    
%     bar(0:1, fcon*100)
    bar(fcon*100, 'stacked')
%     errorbar(0:1, fcon*100, se_fcon*100/2, 'k', 'LineStyle','none')
    ylim([0 105])
    ylabel('% Continuing constriction')
    set(gcf, 'UserData', ud)
end