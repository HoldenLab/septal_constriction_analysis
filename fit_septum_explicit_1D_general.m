% Author: Kevin Whitley
% Date created: 190723

% This function fits an image of a bacterial septum to an explicit model
% using a minimization routine.

function [improf, orthprof, fitvals, fitvals_ax] = fit_septum_explicit_1D_general(imstack, plot_im, param, xy0, olaystack_cropped)

blobedge = 520; % [nm] for cutting out assuming a square. This value seems to work well. Updated 200217
blobarea_pix = blobedge^2 / param.pixSz^2; % [pix^2]

r0 = 2000 / 2 / param.pixSz; % [pix] half of line length for profile. updated 200217
ybox = floor(500 / 2 / param.pixSz); % [pix] width of septum for line profile
xbox = floor(1000 / 2 / param.pixSz); % [pix] width of septum for axial line profile. also the name of a popular game console.
r0r = floor(1500 / 2 / param.pixSz); % [pix] half of line length for profile.
r0ax = floor(2000 / 2 / param.pixSz); % [pix] half of line length for profile.

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

superGaussian = @(a,x) a(4) * exp(-((x-a(1)).^2 / (2*a(2)^2)).^ a(3)) + a(5);
%         superGaussian = @(a,x)  exp(-((x-a(1)).^2 / (2*a(2)^2)).^ a(3));

improf=[]; orthprof=[]; fitvals=[]; fitvals_ax=[];
for ii = 1:size(imstack,3)
    
    frame = imstack(:,:,ii);
    
    if param.plot_gauss
        hold(h_gauss, 'off')
    end
    if param.plot_explicit
        hold(h_exp, 'off')
    end
    
    %% Process image - find septal axis
    
    % WITH SEGMENTED IMAGE
    if ~isempty(olaystack_cropped)
        % lets calculate the imrot angles and centroids and eliminate outliers/ failures
        %         t = 1:size(imstack,3);
        %         theta = 0.*t;
        bactOutline0 = olaystack_cropped(:,:,ii);
        %isolate only the bact blob marked by current centroid
        bactOutline = bwselect(bactOutline0,xy0(ii,1),xy0(ii,2),4);
        %TODO : what if we get more than 1 blob?
        props = regionprops(bactOutline, 'Orientation');
        if ~isempty(props)
            theta = (-props.Orientation +90) * pi/180; % minus sign because it needs to be reflected about y axis.
        else
            %             theta0 = NaN;
            improf(ii,:) = ones(1,size(improf,2))*NaN;
            orthprof(ii,:) = ones(size(orthprof,2),1)*NaN;
            fitvals(ii,:) = [NaN NaN NaN];
            fitvals_ax(ii,:) = [NaN NaN];
%             fitvals_ax(ii,:) = [NaN NaN NaN];
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
            improf(ii,:) = ones(1,size(improf,2))*NaN;
            orthprof(ii,:) = ones(size(orthprof,2),1)*NaN;
            fitvals(ii,:) = [NaN NaN NaN];
            fitvals_ax(ii,:) = [NaN NaN];
%             fitvals_ax(ii,:) = [NaN NaN NaN];
            continue
        end
    end
    
    %% Get intensity profile across septum
    
    % rotate image by theta (x = septal axis, y = cell axis). uses
    % nearest-neighbor interpolation rather than bilinear because black,
    % cropped edges will be artificially filled in with bilinear
    % interpolation (fake signal).
%     theta = theta + pi/2; % HACK!!!!!
    rotim = imrotate(frame, theta*180/pi);
    
    % find new center of septum (moves after imrotate operation)
    Rz = rotz(-theta*180/pi); % don't know why it needs to be -theta, but it does
    vec1 = [xy0(1)-(size(frame,1)+1)/2; xy0(2)-(size(frame,2)+1)/2; 0]; % vector from center of septum to center of image
    vec1rot = Rz*vec1; % rotated vector
    xy0r = [vec1rot(1)+(size(rotim,1)+1)/2; vec1rot(2)+(size(rotim,2)+1)/2; 0]; % new center of septum (image rotated, size changed)
    xy0r = round(xy0r); % image is discrete pixels - need integer values

    % crop image to a thin rectangle along septal axis
    yrange_sepim = max([1 xy0r(2)-ybox]):min([size(rotim,1) xy0r(2)+ybox]);
    xrange_sepim = max([1 xy0r(1)-r0r]):min([size(rotim,1) xy0r(1)+r0r]);
    sepim = rotim(yrange_sepim, :); % crop - only image of septum +/- a few pixels
%     sepim = rotim(yrange_sepim, xrange_sepim); % crop - only image of septum +/- a few pixels
    sepim(sepim==0) = NaN; % zeros are from rotation, and are artificial
    
    % take average of all line profiles across septal axis
    ip = nanmean(sepim,1);
    
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
    axim = rotim(:, xrange_axim); % crop - only image of septum +/- a few pixels
    axim(axim==0) = NaN; % zeros are from rotation, and are artificial
    
    ip_orth = nanmean(axim,2)';
    
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
%         fitvals_ax(ii,:) = [NaN NaN NaN];
        continue
    end
    
    imp = improf(ii,~isnan(improf(ii,:)));
%     imp_proc = (imp-min(imp));
%     imp_proc = imp_proc/max(imp_proc);
    imp_proc = imp;
    [pks, locs] = findpeaks(imp_proc);
    pkmat = [pks' locs'];
    pkmat = sortrows(pkmat);
    if size(pkmat,1)>1 && pkmat(end,1)-pkmat(end-1,1)<0.9
        initwidth = abs((pkmat(end,2)-pkmat(end-1,2))/1.5);
    else
        initwidth = 3;
    end
%     initguess = [length(imp)/2+1 initwidth];
    initguess = [length(imp)/2+1 initwidth min(imp) max(imp)];
    options = optimset('Display', 'off');
%     range = 1:0.1:length(imp);
    range = 0.5:0.1:length(imp)+0.5;
    [fitex, fval] = fminsearch(@(x)septum_1D_obj_amp(imp,x,range,param), initguess, options); % minimize residual sum of squares
%     [fitex, fval] = fminsearch(@(x)septum_1D_obj(imp_proc,x,range,param), initguess, options);
%     [fitex, fval] = fmincon(@(x)septum_1D_obj(imp_proc,x,range,param), initguess, [],[],[],[],[0 0],[Inf Inf]);

    sst_rad = sum((imp - mean(imp)).^2);
    Rsq_rad = 1 - fval/sst_rad;
    
%     fitvals(ii,:) = [fitex(1) fitex(2) fval];
    fitvals(ii,:) = [fitex(1) fitex(2) Rsq_rad];
    
    if param.plot_explicit
%         plot(h_exp, range(1):length(imp_proc), imp_proc)
        plot(h_exp, (range(1):length(imp_proc))-fitex(1), imp_proc)
        hold(h_exp, 'on')

        plot(h_exp, range(1:end-1)-fitex(1), septum_model_1D_amp(fitvals(ii,1),fitvals(ii,2),param.psfFWHM,param.pixSz,range,fitvals(ii,3),fitvals(ii,4)))
%         plot(h_exp, range(1:end-1)-fitex(1), septum_model_1D(fitvals(ii,1),fitvals(ii,2),param.psfFWHM,param.pixSz,range))
%         figure(hRawInt);
%         plot(imp);
        ylabel('raw intensity')
    end
    
    %% Fit orthogonal intensity line profile to generalized (or 'super') Gaussian model
    
    linep_ax = orthprof(ii,~isnan(orthprof(ii,:)));
    if isempty(linep_ax) || length(linep_ax)<5 || ~any(linep_ax~=linep_ax(1))
        FWHM_ax = NaN;
%         se_ax = ones(1,5)*NaN;
    else
        linep_ax_proc = linep_ax - min(linep_ax);
        linep_ax_proc = linep_ax_proc / max(linep_ax_proc);
        h = size(linep_ax_proc, 2);
        x = 1:h;
        
%         if param.plot_gauss
%             plot(h_gauss, linep_ax)
%             hold(h_gauss, 'on')
%         end
        
        initwidth = 5;
%         initguess = [mean(x), initwidth, 1];
        initguess_g = [mean(x) initwidth 1 max(linep_ax) min(linep_ax)];
%         initguess = [mean(x) initwidth min(linep_ax) max(linep_ax)];
%         options = optimoptions('lsqcurvefit', 'Display', 'off');
%         lb = [0 0 0.5];
        lb = [0 0 1 0 0];
        [a, rn, res, ~, ~, ~, J] = lsqcurvefit(superGaussian, initguess_g, x, linep_ax, lb, [], options);
%         range = 0.5:0.1:length(linep_ax)+0.5;
%         [fitax, fval_ax] = fminsearch(@(x)septum_thickness_1D_obj_amp(linep_ax,x,range,param), initguess, options); % minimize residual sum of squares
        sst_ax = sum((linep_ax - mean(linep_ax)).^2);
        Rsq_ax = 1 - rn/sst_ax;
        
%         [~, se_ax] = nlparci2(a, real(res), 'Jacobian', real(J));

%         xC = a(1);
        par_d = a(2);
        par_e = a(3);
        FWHM_ax = 2*sqrt(2)*par_d*(log(2))^(1/(2*par_e));
%         x0ax = xC - FWHM_ax/2;
%         x1 = xC + FWHM_ax/2;
%         widthResult = [xC, x0ax, x1];
%         intensity = sum(linep_ax_proc);
        
        if param.plot_gauss
            plot(h_gauss, (range(1):length(linep_ax)), linep_ax)
            hold(h_gauss, 'on')
            tv = x(1):(x(2)-x(1))/10:x(end);
            plot(h_gauss, tv-fitax(1), superGaussian(a,tv))
%             plot(h_gauss, range(1:end-1)-fitax(1), septum_thickness_model_1D_amp(fitax(1),fitax(2),param.psfFWHM,param.pixSz,range,fitax(3),fitax(4)))
        end

    end
    
%     fitvals_ax(ii,:) = [FWHM_ax se_ax(2)];
    fitvals_ax(ii,:) = [FWHM_ax Rsq_ax];
%     fitvals_ax(ii,:) = [FWHM_ax rn];
%     fitvals_ax(ii,:) = [fitax(1) fitax(2) fval_ax];
%     fitvals_ax(ii,:) = [fitax(1) fitax(2) Rsq_ax];
    
    %% Plot image with septal and axial axes
    
    if plot_im
        hold(h_im, 'off')
        imagesc(h_im, rotim)
        
        hold(h_im, 'on')
        xs = xy0r(1)-r0:0.1:xy0r(1)+r0;
        plot(h_im, xs, ones(1,length(xs))*xy0r(2), 'k', 'linew', 2)
        ys = xy0r(2)-r0:0.1:xy0r(2)+r0;
        plot(h_im, ones(1,length(ys))*xy0r(1), ys, 'r', 'linew', 2)
%         figure(hIm0);
%         imagesc(frame);
%         figure(hIm1);
%         imagesc(seg_f);
%         ii
%         pause
    end
end
