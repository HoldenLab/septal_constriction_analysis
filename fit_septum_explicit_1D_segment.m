% Author: Kevin Whitley
% Date created: 190723

% This function fits an image of a bacterial septum to an explicit model
% using a minimization routine.

function [improf, orthprof, fitvals, fitvals_ax] = fit_septum_explicit_1D_segment(imstack, plot_im, param, xy0, olaystack_cropped)

plot_gauss = param.plot_gauss;
plot_explicit = param.plot_explicit;

blobedge = 520; % [nm] for cutting out assuming a square. This value seems to work well. Updated 200217
blobarea_pix = blobedge^2 / param.pixSz^2; % [pix^2]

%r0 = 2000 / 2 / param.pixSz; % [pix] half of line length for profile. updated 200217
r0r = floor(1500 / 2 / param.pixSz); % [pix] half of line length for profile.
r0ax = floor(2000 / 2 / param.pixSz); % [pix] half of line length for profile.
ybox = floor(500 / 2 / param.pixSz); % [pix] width of septum for line profile
xbox = floor(1000 / 2 / param.pixSz); % [pix] width of septum for axial line profile. also the name of a popular game console.

if plot_im
    figure('Position',[1100 450 400 300])
    h_im = gca;
end
if param.plot_gauss
    figure('Position',[700 450 400 300])
    h_gauss = gca;
end
if param.plot_explicit
    figure('Position',[700 50 400 300])
    h_exp = gca;
end

%before we do the fit
% lets calculate the imrot angles and centroids and eliminate outliers/ failures
t = 1:size(imstack,3);
theta = 0.*t;
for ii = 1:size(imstack,3)
    im = imstack(:,:,ii);
    bactOutline0 = olaystack_cropped(:,:,ii);
    %isolate only the bact blob marked by current centroid
    bactOutline = bwselect(bactOutline0,xy0(ii,1),xy0(ii,2),4);
    %TODO : what if we get more than 1 blob? 
    props = regionprops(bactOutline, 'Orientation');
    if ~isempty(props)
        theta0(ii)= (-props.Orientation +90) * pi/180; % minus sign because it needs to be reflected about y axis.
    else
        theta0(ii) = NaN;
    end
end 
%we could fit a line to remove outliers from theta
%use a robust fit
% brob=robustfit(t,theta0);
% theta=brob(1)+brob(2)*t;
theta = theta0;
%%<DEBUG
%figure;
%hold all;
%plot(t,theta0);
%plot(t,theta);
%ylabel('theta');
%pause
%%DEBUG>

superGaussian = @(a,x)  a(4)*exp(-((x-a(1)).^2 / (2*a(2)^2)).^ a(3))+a(5);

improf=[]; orthprof=[]; fitvals=[]; fitvals_ax=[];
for ii = 1:size(imstack,3)
    
    frame = imstack(:,:,ii);
    
    if param.plot_gauss
        hold(h_gauss, 'off')
    end
    if param.plot_explicit
        hold(h_exp, 'off')
    end
    
    % multilevel thresh
    thr = multithresh(frame, 3);
    seg_f = imquantize(frame, thr);
    
    bwframe = seg_f > length(thr); % top level
    
    if ~isfield(param,'ZGFP') || param.ZGFP == 0
        cc = bwconncomp(bwframe);
        if cc.NumObjects < 2
            bwframe2 = seg_f > length(thr)-1;
            bwframe = bwframe2;
        end
    end
    
    % remove all but one or two blobs from image
    BW = bwareafilt(bwframe, 2, 'largest');
    bwframe(BW==0) = 0;

    %% Process image - try to find two blobs
    
    if isnan(theta(ii))
        improf(ii,:) = ones(1,size(improf,2))*NaN;
        orthprof(ii,:) = ones(size(orthprof,2),1)*NaN;
        fitvals(ii,:) = [NaN NaN NaN NaN];
        fitvals_ax(ii,:) = [NaN NaN];
        continue
    end
    
%     props2 = regionprops(bwframe, 'Centroid', 'Area', 'Orientation', 'MajorAxisLength', 'Area');
%     cents = cat(1, props2.Centroid); % centroids in image
    
%     if size(cents,1) == 2 % nascent ring
%         % exception: if the area of either blob is too big, it's basically
%         % never a nascent ring.
%         if props2(1).Area>blobarea_pix || props2(2).Area>blobarea_pix
%             
%             % mature ring. added 200217
%             if props2(1).Area > props2(2).Area
%                 ind = 1;
%             else
%                 ind = 2;
%             end
%             xy0 = cents(ind,:);
%             %         theta = -props(ind).Orientation * pi/180; % minus sign because it needs to be reflected about y axis.
%         else
%             width = (cents(2,:) - cents(1,:)) / 2; % distance from centroid to centroid
%             xy0 = cents(2,:) - width; % position of center between centroids
%             %         theta = atan((cents(2,2)-cents(1,2))./(cents(2,1)-cents(1,1))); % orientation of septum
%         end
%     elseif size(cents,1) == 1 % mature ring
%         xy0 = cents;
%         %     theta = -props.Orientation * pi/180; % minus sign because it needs to be reflected about y axis.
%     else % probably no centroid because image is utter crap
%         improf(ii,:) = ones(1,size(improf,2))*NaN;
%         orthprof(ii,:) = ones(size(orthprof,2),1)*NaN;
%         fitvals(ii,:) = [NaN NaN NaN NaN];
%         fitvals_ax(ii,:) = [NaN NaN];
%         continue
%     end
    %% Get intensity profile across septum
    
    % rotate image by theta (x = septal axis, y = cell axis)
    rotim = imrotate(frame, theta(ii)*180/pi);
    if plot_im
        hold(h_im, 'off')
        imagesc(h_im, rotim)
    end
    
    Rz = rotz(-theta(ii)*180/pi); % don't know why it needs to be -theta, but it does
%     vec1 = [xy0(ii,1)-(size(frame,1)+1)/2; xy0(ii,2)-(size(frame,2)+1)/2; 0]; % vector from center of septum to center of image
    vec1 = [xy0(1)-(size(frame,1)+1)/2; xy0(2)-(size(frame,2)+1)/2; 0]; % vector from center of septum to center of image
    vec1rot = Rz*vec1; % rotated vector
    xy0r = [vec1rot(1)+(size(rotim,1)+1)/2; vec1rot(2)+(size(rotim,2)+1)/2; 0]; % new center of septum (image rotated, size changed)
    xy0r = round(xy0r); % image is discrete pixels - need integer values
    
    yrange_sepim = max([1 xy0r(2)-ybox]):min([size(rotim,1) xy0r(2)+ybox]);
    xrange_sepim = max([1 xy0r(1)-r0r]):min([size(rotim,1) xy0r(1)+r0r]);
    sepim = rotim(yrange_sepim, xrange_sepim); % crop - only image of septum +/- a few pixels
    sepim(sepim==0) = NaN; % zeros are from rotation, and are artificial
    
    %ip = nanmean(sepim,1);
    ip = nansum(sepim,1);
    
    % pad ip or improf with NaNs if necessary
    if ~isempty(improf) && length(ip)>length(improf(end,:))
        ldiff = length(ip) - length(improf(end,:));
        newcol = ones(size(improf,1),ldiff)*NaN;
        improf = [improf newcol];
    elseif ~isempty(improf) && length(ip)<length(improf(end,:))
        ip = [ip ones(1,length(improf(end,:))-length(ip))*NaN];
    end
    
    if ii>1 && isempty(improf)
        improf = ones(ii-1,length(ip))*NaN;
    end
    
    improf(ii,:) = ip;
    
    %% Get intensity profile orthogonal to septum
    
    xrange_axim = max([1 xy0r(1)-xbox]):min([size(rotim,2) xy0r(1)+xbox]);
    yrange_axim = max([1 xy0r(2)-r0ax]):min([size(rotim,2) xy0r(2)+r0ax]);
    axim = rotim(yrange_axim, xrange_axim); % crop - only image of septum +/- a few pixels
    axim(axim==0) = NaN; % zeros are from rotation, and are artificial
    
    %ip_orth = nanmean(axim,2)';
    ip_orth = nansum(axim,2)';
    
    % pad ip_orth or orthprof with NaNs if necessary
    if ~isempty(orthprof) && size(ip_orth,2)>size(orthprof,2)
        ldiff = size(ip_orth,2) - size(orthprof,2);
        newcol = ones(size(orthprof,1),ldiff)*NaN;
        orthprof = [orthprof newcol];
    elseif ~isempty(orthprof) && size(ip_orth,2)<size(orthprof,2)
        ip_orth = [ip_orth ones(size(ip_orth,1),size(orthprof,2)-size(ip_orth,2))*NaN];
    end
    
    if ii>1 && isempty(orthprof)
        orthprof = ones(ii-1,size(ip_orth,2))*NaN;
    end
    
    orthprof(ii,:) = ip_orth;

    %% Fit septal intensity line profile to explicit 'tilted circle' model
    
    % line profile is too short for fitting - ignore
    if length(ip(~isnan(ip)))<3
        fitvals(ii,:) = [NaN NaN NaN];
        fitvals_ax(ii,:) = [NaN NaN];
        continue
    end
    
    imp = improf(ii,~isnan(improf(ii,:)));
    imp_proc = (imp-min(imp));
    imp_proc = imp_proc/max(imp_proc);
    [pks, locs] = findpeaks(imp_proc);
    pkmat = [pks' locs'];
    pkmat = sortrows(pkmat);
    if size(pkmat,1)>1 && pkmat(end,1)-pkmat(end-1,1)<0.9
        initwidth = abs((pkmat(end,2)-pkmat(end-1,2))/1.5);
    else
        initwidth = 3;
    end
    initguess = [length(imp)/2+1 initwidth max(imp_proc)];
    options = optimset('Display', 'off');
    range = 1:0.1:length(imp);
    [fitex, fval] = fminsearch(@(x)septum_1D_objZ(imp_proc,x,range,param), initguess, options);
    
    fitvals(ii,:) = [fitex fval];
     
     if param.plot_explicit
         halfrange = range - (range(2)-range(1))/2;
         plot(h_exp, imp_proc)
         hold(h_exp, 'on')
         plot(h_exp, halfrange, septum_model_1DZ(fitvals(ii,1),fitvals(ii,2),250,65,range, fitvals(ii,3)));
%          figure(hRawInt);
%          plot(imp);
%          ylabel('raw intensity')
     end

    %%alternative supergaussian septal fit
    %h = length(imp_proc);
    %x = 1:h;
    %initguess = [mean(x), initwidth, 1,1,0];
    %options = optimoptions('lsqcurvefit', 'Display', 'off');
    %lb = [0 0 0.5 0 0];
    %[a, rn, res, ~, ~, ~, J] = lsqcurvefit(superGaussian, initguess, x, imp_proc, lb, [], options);
    %[~, se_rad] = nlparci2(a, real(res), 'Jacobian', real(J));

    %xC = a(1);
    %par_d = a(2);
    %par_e = a(3);
    %FWHM_rad = 2*sqrt(2)*par_d*(log(2))^(1/(2*par_e));
    %
    %if param.plot_explicit
    %    tv = x(1):(x(2)-x(1))/10:x(end);
    %    hold(h_exp, 'off')
    %    plot(h_exp, imp_proc)
    %    hold(h_exp, 'on')
    %    plot(h_exp, tv, superGaussian(a,tv))
    %end
    %fitvals(ii,:) = [FWHM_rad se_rad(2)];

    %% Fit orthogonal intensity line profile to generalized (or 'super') Gaussian model
    
    linep_ax = orthprof(ii,~isnan(orthprof(ii,:)));
    if isempty(linep_ax) || length(linep_ax)<5
        FWHM_ax = NaN;
        se_ax = ones(1,5)*NaN;
    else
        linep_ax_proc = linep_ax - min(linep_ax);
        linep_ax_proc = linep_ax_proc / max(linep_ax_proc);
        h = size(linep_ax_proc, 2);
        x = 1:h;
        
        if param.plot_gauss
            plot(h_gauss, linep_ax_proc)
            hold(h_gauss, 'on')
        end

        initwidth = 5;
        initguess = [mean(x), initwidth, 1,1,0];
        options = optimoptions('lsqcurvefit', 'Display', 'off');
        lb = [0 0 0.5 0 0];
        [a, rn, res, ~, ~, ~, J] = lsqcurvefit(superGaussian, initguess, x, linep_ax_proc, lb, [], options);
        sst = sum((linep_ax_proc - mean(linep_ax_proc)).^2);
        Rsq = 1 - rn/sst;
        
        [~, se_ax] = nlparci2(a, real(res), 'Jacobian', real(J));

        xC = a(1);
        par_d = a(2);
        par_e = a(3);
        FWHM_ax = 2*sqrt(2)*par_d*(log(2))^(1/(2*par_e));
        x0ax = xC - FWHM_ax/2;
        x1 = xC + FWHM_ax/2;
        widthResult = [xC, x0ax, x1];
        intensity = sum(linep_ax_proc);
        
        if param.plot_gauss
            tv = x(1):(x(2)-x(1))/10:x(end);
            plot(h_gauss, tv, superGaussian(a,tv))
        end

    end
    
    fitvals_ax(ii,:) = [FWHM_ax se_ax(2)];
    
    %% Plot image with septal and axial axes
    
    if plot_im
%         hold(h_im, 'on')
        
        hold(h_im, 'on')
        xs = xy0r(1)-r0r:0.1:xy0r(1)+r0r;
        plot(h_im, xs, ones(1,length(xs))*xy0r(2), 'k', 'linew', 2)
        ys = xy0r(2)-r0ax:0.1:xy0r(2)+r0ax;
        plot(h_im, ones(1,length(ys))*xy0r(1), ys, 'r', 'linew', 2)
%         ii
%         pause
    end
end
%-----------------------------------------------------------
% Author: Kevin Whitley
% Date created: 200128

% This is an objective function to minimize for fitting
% horizontally-oriented bacterial septal profiles.

function obj = septum_1D_objZ(prof_im, guess, range, param)

if nargin == 0
    guess = [10 15 1 0];
    range = 1:length(prof_im);
end

mu = guess(1);
R = guess(2);
% ring_grad = guess(3);
amp = guess(3);
psfFWHM = param.psfFWHM;
pixSz = param.pixSz;

prof_model = septum_model_1DZ(mu, R, psfFWHM, pixSz, range,amp);

sp = range(2)-range(1); % spacing
prof_model = downsample(prof_model, uint8(1/sp)); % same number of points as image profile

diff = prof_im - prof_model;

obj = sum(sum(diff.^2));
%----------------------------------------------------------------
% Author: Kevin Whitley
% Date created: 200128

% This function produces a line profile of a septum using an explicit
% 'tilted circle' model.

function improf = septum_model_1DZ(X0, R, psfFWHM, pixSz, X,amp)

if nargin==0
    X0 = 13.1;
    R = 6; % [pix]
    ring_grad = 0;
    psfFWHM = 250;
    pixSz = 65;
    sp = .01; % spacing needs to be low enough to prevent problems with discretization (probably <0.5 or so)
    X = 1:sp:25;
    amp=10;
end

sigma = psfFWHM/2.35/pixSz;

x1 = X(1:end-1) - X0;
x2 = X(2:end) - X0;
% arcl_i = 2*R * real(asin(x1/R)) + (1:length(x1))*ring_grad; % arc length out to axi
arcl_i = 2*R * real(asin(x1/R)); % arc length out to axi
arcl_ip1 = 2*R * real(asin(x2/R)); % arc length out to ax(i+1)

ring_prof = abs(arcl_ip1 - arcl_i); % Delta(arc length)

% ring_prof(abs(x1)>R) = 0;

gpdf_fcn = @(a,x) 1/(a(2)*sqrt(2*pi)) * exp(-(x-a(1)).^2 ./ (2*a(2)^2));
gauss = gpdf_fcn([0 sigma], x1);

ring_prof = conv(ring_prof, gauss); % convolve with Gaussian psf
ring_prof = ring_prof(floor(length(gauss)/2):floor((end-length(gauss)/2))+1); % cut ends (convolution increases length of vector)

ring_prof_scale = ring_prof/max(ring_prof(:));

improf = ring_prof_scale*amp;

if nargin==0
    figure
    hold on
    halfX = X - (X(2)-X(1))/2;
    plot(halfX, improf)
end

