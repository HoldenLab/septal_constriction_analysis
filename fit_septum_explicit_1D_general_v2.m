% Author: Kevin Whitley
% Date created: 190723

% This function fits an image of a bacterial septum to an explicit model
% using a minimization routine.

% v2 (201204):
%   - initialized variables better
%   - cleaned up a bunch of crap
%   - made objective and model functions private

function [improfs, orthprof, fitvals, fitvals_ax] = fit_septum_explicit_1D_general_v2(imstack, param, xy0, olaystack_cropped)

blobedge = 520; % [nm] for cutting out assuming a square. This value seems to work well. Updated 200217
blobarea_pix = blobedge^2 / param.pixSz^2; % [pix^2]

r0 = 2000 / 2 / param.pixSz; % [pix] half of line length for profile. updated 200217
ybox = floor(500 / 2 / param.pixSz); % [pix] width of septum for line profile
xbox = floor(1000 / 2 / param.pixSz); % [pix] width of septum for axial line profile. also the name of a popular game console.
r0ax = floor(1200 / 2 / param.pixSz); % [pix] half of line length for profile.

if param.plot_im_fits
    figure('Position',[700 100 800 600])
    h_im = subplot(221);
    h_rotim = subplot(222);
    h_rad = subplot(223);
    h_ax = subplot(224);
end

superGaussian = @(a,x) a(4) * exp(-((x-a(1)).^2 / (2*a(2)^2)).^ a(3)) + a(5); % background and amplitude float

% initialize output variables with NaNs
improfs = NaN(size(imstack,3), ceil(size(imstack,1)*sqrt(2)));
orthprof = improfs;
fitvals = NaN(size(imstack,3),3);
fitvals_ax = NaN(size(imstack,3),2);

for ii = 1:size(imstack,3)
    %% Process image - find septal axis
    
    frame = imstack(:,:,ii);
    
    % WITH SEGMENTED IMAGE
    if ~isempty(olaystack_cropped)

        % isolate only the bact blob marked by current centroid
        bactOutline0 = olaystack_cropped(:,:,ii);
        bactOutline = bwselect(bactOutline0,xy0(ii,1),xy0(ii,2),4);
        
        % TODO : what if we get more than 1 blob?
        
        props = regionprops(bactOutline, 'Orientation');
        if ~isempty(props)
            theta = (-props.Orientation+90) * pi/180; % minus sign because it needs to be reflected about y axis.
        else
            continue
        end
    else
        
    % WITHOUT SEGMENTED IMAGE
     
    % FIND TWO BLOBS IF POSSIBLE
        
        % multilevel thresh
        thr = multithresh(frame, 3);
        seg_f = imquantize(frame, thr);
        
        bwframe = seg_f > length(thr); % top level
        
        % if bwframe only has one object, drop down one level and try again
        % (helps with asymmetric nascent septa). Seems mostly to help with
        % 2B-mNG strain, not so much with Z-GFP strain.
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
        
        % FIND SEPTAL AXIS
        
        props = regionprops(bwframe, 'Centroid', 'Area', 'Orientation', 'MajorAxisLength', 'Area');
        cents = cat(1, props.Centroid); % centroids in image
        
        if size(cents,1) == 2 % nascent ring
            
            % exception: if the area of either blob is too big, it's basically
            % never a nascent ring.
            if props(1).Area>blobarea_pix || props(2).Area>blobarea_pix
                
                % mature ring. added 200217
                if props(1).Area > props(2).Area
                    ind = 1;
                else
                    ind = 2;
                end
                xy0 = cents(ind,:);
                theta = -props(ind).Orientation * pi/180; % minus sign because it needs to be reflected about y axis.
            else
                width = (cents(2,:) - cents(1,:)) / 2; % distance from centroid to centroid
                xy0 = cents(2,:) - width; % position of center between centroids
                theta = atan((cents(2,2)-cents(1,2))./(cents(2,1)-cents(1,1))); % orientation of septum
            end
            
        elseif size(cents,1) == 1 % mature ring
            
            xy0 = cents;
            theta = -props.Orientation * pi/180; % minus sign because it needs to be reflected about y axis.
            
        else % probably no centroid because image is utter crap
            continue
        end
    end
    
    %% Get intensity profiles
    
    % ROTATE IMAGE
    
    % rotate image by theta (x = septal axis, y = cell axis). uses nearest-neighbor interpolation rather than bilinear because black,
    % cropped edges will be artificially filled in with bilinear interpolation (fake signal).
    rotim = imrotate(frame, theta*180/pi);
    
    % find new center of septum (moves after imrotate operation)
    Rz = rotz(-theta*180/pi); % don't know why it needs to be -theta, but it does
    vec1 = [xy0(1)-(size(frame,1)+1)/2; xy0(2)-(size(frame,2)+1)/2; 0]; % vector from center of septum to center of image
    vec1rot = Rz*vec1; % rotated vector
    xy0r = [vec1rot(1)+(size(rotim,1)+1)/2; vec1rot(2)+(size(rotim,2)+1)/2; 0]; % new center of septum (image rotated, size changed)
    xy0r = round(xy0r); % image is discrete pixels - need integer values

    % INTENSITY PROFILE ACROSS SEPTUM
    
    % crop image to a thin rectangle along septal axis
    yrange_sepim = max([1 xy0r(2)-ybox]):min([size(rotim,1) xy0r(2)+ybox]);
    sepim = rotim(yrange_sepim, :); % crop - only image of septum +/- a few pixels
    sepim(sepim==0) = NaN; % zeros are from rotation, and are artificial
    
    % take average of all line profiles across septal axis for this frame
    ip = nanmean(sepim,1);

    % add single-frame line profile to matrix of line profiles
    improfs(ii,1:length(ip)) = ip;
    
    % INTENSITY PROFILE ORTHOGONAL TO SEPTUM
    
    xrange_axim = max([1 xy0r(1)-xbox]):min([size(rotim,2) xy0r(1)+xbox]);
    yrange_axim = max([1 xy0r(2)-r0ax]):min([size(rotim,2) xy0r(2)+r0ax]);
    axim = rotim(yrange_axim, xrange_axim); % crop - only image of septum +/- a few pixels. Updated 201029 to have yrange_axim
    axim(axim==0) = NaN; % zeros are from rotation, and are artificial
    
    % take average of all line profiles orthogonal to septal axis for this frame
    ip_ax = nanmean(axim,2)';
    
    % add single-frame line profile to matrix of line profiles
    orthprof(ii,1:length(ip_ax)) = ip_ax;

    %% Fit line profiles to appropriate models
    
    % FIT LINE PROFILE ALONG SEPTUM
    
    ip = ip(~isnan(ip));
    
    if length(ip)>=3

        % initial guess parameters
        [pks, locs] = findpeaks(ip);
        pkmat = [pks' locs'];
        pkmat = sortrows(pkmat);
        if size(pkmat,1)>1 && pkmat(end-1,1)/pkmat(end,1)>=0.9 % two peaks with similar heights
            initwidth = abs((pkmat(end,2)-pkmat(end-1,2))/1.5);
        else % only one peak, or very asymmetric
            initwidth = 5;
        end
        initguess = [length(ip)/2+1 initwidth min(ip) max(ip)];
        
        % fit line profile to tilted circle model
        %     range = 1:0.1:length(imp);
        x_rad = 0.5:0.1:length(ip)+0.5;
        options = optimset('Display', 'off');
        [fitex, fval] = fminsearch(@(x)septum_1D_obj_amp_private(ip,x,x_rad,param), initguess, options); % minimize residual sum of squares
        
        % calculate fitting error
        sst_rad = sum((ip - mean(ip)).^2);
        Rsq_rad = 1 - fval/sst_rad;
        
        %     fitvals(ii,:) = [fitex(1) fitex(2) fval]; % changed back 201028
        fitvals(ii,:) = [fitex(1) fitex(2) Rsq_rad];
        
    end
    
    % FIT LINE PROFILE ORTHOGONAL TO SEPTUM
    
    ip_ax = ip_ax(~isnan(ip_ax));
    
    if length(ip_ax)>=3 && any(ip_ax~=ip_ax(1)) % need enough points, and not a flat line
 
        % initial guess parameters
        initwidth = 5;
        initguess_ax = [length(ip_ax)/2+1 initwidth 1 max(ip_ax) min(ip_ax)];
        
        % fit line profile to supergaussian function, calculate FWHM
        x_ax = 1:length(ip_ax);
        lb = [0 0 1 0 0];
        [fitax, rn] = lsqcurvefit(superGaussian, initguess_ax, x_ax, ip_ax, lb, [], options);
        
        sigma = fitax(2);
        superexp = fitax(3);
        FWHM_ax = 2*sqrt(2)*sigma*(log(2))^(1/(2*superexp));
        
        % calculate fitting error
        sst_ax = sum((ip_ax - mean(ip_ax)).^2);
        Rsq_ax = 1 - rn/sst_ax;
        
        fitvals_ax(ii,:) = [FWHM_ax Rsq_ax];
        
    end

%% Plot images, line profiles, and fits
    
    if param.plot_im_fits
        hold(h_im, 'off')
        imagesc(h_im, frame)
        title(h_im, 'Raw image')
        
        hold(h_rotim, 'off')
        imagesc(h_rotim, rotim)
        hold(h_rotim, 'on')
        
        xs = xy0r(1)-r0:0.1:xy0r(1)+r0;
        plot(h_rotim, xs, ones(1,length(xs))*xy0r(2), 'k', 'linew', 2)
        ys = xy0r(2)-r0:0.1:xy0r(2)+r0;
        plot(h_rotim, ones(1,length(ys))*xy0r(1), ys, 'r', 'linew', 2)
        title(h_rotim, 'Rotated image')

        hold(h_rad, 'off')
        plot(h_rad, (x_rad(1):length(ip))-fitex(1), ip, 'Color', 'k', 'linew', 2)
        hold(h_rad, 'on')

        plot(h_rad, x_rad(1:end-1)-fitex(1), septum_model_1D_amp_private(fitex(1),fitex(2),param.psfFWHM,param.pixSz,x_rad,fitex(3),fitex(4)), 'b')
        title(h_rad, 'Line profile along septum')
        xlabel(h_rad, 'Distance along septum (pix)')
        ylabel(h_rad, 'Raw intensity')

        hold(h_ax, 'off')
        plot(h_ax, (x_ax(1):length(ip_ax)), ip_ax, 'Color', 'r', 'linew', 2)
        hold(h_ax, 'on')
        
        tv = x_ax(1):(x_ax(2)-x_ax(1))/10:x_ax(end);
        plot(h_ax, tv, superGaussian(fitax,tv), 'b')
        title(h_ax, 'Line profile orthogonal to septum')
        xlabel(h_ax, 'Distance orthogonal to septum (pix)')
        ylabel(h_ax, 'Raw intensity')
    end
end

function obj = septum_1D_obj_amp_private(prof_im, guess, range, param)
% This is an objective function to minimize for fitting
% horizontally-oriented bacterial septal profiles.

mu = guess(1);
R = guess(2);
bg = guess(3);
amp = guess(4);
psfFWHM = param.psfFWHM;
pixSz = param.pixSz;

prof_model = septum_model_1D_amp_private(mu, R, psfFWHM, pixSz, range, bg, amp);

sp = range(2)-range(1); % spacing
prof_model = downsample(prof_model, uint8(1/sp)); % same number of points as image profile

diff = prof_im - prof_model;

obj = sum(diff.^2);


function improf = septum_model_1D_amp_private(X0, R, psfFWHM, pixSz, X, bg, amp)
% This function produces a line profile of a septum using an explicit
% 'tilted circle' model. Discrete version of model.

sigma = psfFWHM/2.35/pixSz;

x1 = X(1:end-1) - X0;
x2 = X(2:end) - X0;
arcl_i = 2*R * real(asin(x1/R)); % arc length out to axi
arcl_ip1 = 2*R * real(asin(x2/R)); % arc length out to ax(i+1)

ring_prof = abs(arcl_ip1 - arcl_i); % Delta(arc length)

gpdf_fcn = @(a,x) 1/(a(2)*sqrt(2*pi)) * exp(-(x-a(1)).^2 ./ (2*a(2)^2));
gauss = gpdf_fcn([0 sigma], x1);

ring_prof = conv(ring_prof, gauss, 'same'); % convolve with Gaussian psf

ring_prof_scale = ring_prof/max(ring_prof(:));

improf = ring_prof_scale * amp + bg;
