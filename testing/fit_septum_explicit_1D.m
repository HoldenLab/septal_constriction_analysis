% Author: Kevin Whitley
% Date created: 190723

% This function fits an image of a bacterial septum to an explicit model
% using a minimization routine.

function [improf, orthprof, fitvals, fitvals_ax] = fit_septum_explicit_1D(imstack, plot_im, param)

plot_gauss = param.plot_gauss;
plot_explicit = param.plot_explicit;

blobedge = 520; % [nm] for cutting out assuming a square. This value seems to work well. Updated 200217
blobarea_pix = blobedge^2 / param.pixSz^2; % [pix^2]

r0 = 2000 / 2 / param.pixSz; % [pix] half of line length for profile. updated 200217

if plot_im
    figure
    h_im = gca;
end
if param.plot_gauss
    figure
    h_gauss = gca;
end
if plot_explicit
    figure('Position',[100 600 500 400])
    h_exp = gca;
end

improf=[]; orthprof=[]; fitvals=[]; fitvals_ax=[];
for ii = 1:size(imstack,3)
    
    frame = imstack(:,:,ii);
    
    if plot_im
        hold(h_im, 'off')
        imagesc(h_im, frame)
    end
    if param.plot_gauss
        hold(h_gauss, 'off')
    end
    if plot_explicit
        hold(h_exp, 'off')
    end
    
    %% Process image - try to find two blobs

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

    %% Find septal axis

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
        continue
    end
    
    %% Get intensity profile across septum

    % endpoints for line going through septal axis
    x0 = [xy0(1)+r0*cos(theta) xy0(1)-r0*cos(theta)];
    y0 = [xy0(2)+r0*sin(theta) xy0(2)-r0*sin(theta)];
    
    ip = improfile(frame, x0, y0)';
    
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
    
    % rotation matrix about z axis
    Rz = rotz(theta*180/pi);
    xy_minv = [min(x0); y0(2); 0]; % starting xy position as 3D vector
    
    orth_space = 0.1; % spacing for each orthogonal line (measured along septal axis)
    ustep = [orth_space; 0; 0]; % unit step size along septal axis
    
    ip_orth = [];
    for jj = 1:r0/orth_space*2
        
        % new 'central' coordinates from rotation matrix (sliding along
        % septal axis)
        new_xy0 = xy_minv + Rz*(ustep.*jj);
        
        new_x0 = new_xy0(1);
        new_y0 = new_xy0(2);
        
        % endpoints for line parallel to cell axis
        x_orth = [new_x0-r0*sin(theta) new_x0+r0*sin(theta)];
        y_orth = [new_y0+r0*cos(theta) new_y0-r0*cos(theta)];

        % line profile along cell axis with shifted 'central' coordinates
        ip_orth(jj,:) = improfile(frame, x_orth, y_orth)';

        if plot_im && 0 % for debugging
            hold(h_im, 'on')
            slp_o1 = (y_orth(2)-y_orth(1)) / (x_orth(2)-x_orth(1));
            yint_o1 = y_orth(2) - (slp_o1*x_orth(2));
            tv_o1 = min([x_orth(1) x_orth(2)]):max([x_orth(1) x_orth(2)]);
            plot(h_im, tv_o1, tv_o1*slp_o1+yint_o1, 'r', 'linew', 3)
            
            plot(h_im, new_xy0(1), new_xy0(2), 'ow');
        end
        if plot_gauss && 0
            plot(h_gauss, ip_orth(jj,:))
        end
    end
    
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
    
%     ip_orth(:,any(isnan(ip_orth))) = NaN; % if any nans are present, sum will not be accurate.
    orthprof(ii,:) = nanmean(ip_orth,1)*orth_space; % sum will be inaccurate due to all the nans present. need to normalize by number of actual data points.
%     orthprof(ii,:) = nansum(ip_orth,1)*orth_space; % integral rather than sum. will not give accurate intensity measurement!

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
    initguess = [length(imp)/2+1 initwidth];
    options = optimset('Display', 'off');
%     range = 1:length(imp);
    range = 1:0.1:length(imp);
    [fitex, fval] = fminsearch(@(x)septum_1D_obj(imp_proc,x,range,param), initguess, options);
    
    fitvals(ii,:) = [fitex fval];
    
    if param.plot_explicit
        halfrange = range - (range(2)-range(1))/2;
        plot(h_exp, imp_proc)
        hold(h_exp, 'on')
        plot(h_exp, halfrange, septum_model_1D(fitvals(ii,1),fitvals(ii,2),250,65,range))
    end
    
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
        
%         superGaussian = @(a,x)  a(1) + a(2)*exp(-((x-a(3)).^2 / (2*a(4)^2)).^ a(5));
%         superGaussian = @(a,x)  a(1)*exp(-((x-a(2)).^2 / (2*a(3)^2)).^ a(4));
        superGaussian = @(a,x)  exp(-((x-a(1)).^2 / (2*a(2)^2)).^ a(3));
        initwidth = 5;
%         initguess = [min(linep_ax_proc), max(linep_ax_proc), mean(x), initwidth, 1];
%         initguess = [max(linep_ax_proc), mean(x), initwidth, 1];
        initguess = [mean(x), initwidth, 1];
        options = optimoptions('lsqcurvefit', 'Display', 'off');
        lb = [0 0 0 0 0.5];
        [a, rn, res, ~, ~, ~, J] = lsqcurvefit(superGaussian, initguess, x, linep_ax_proc, lb, [], options);
        sst = sum((linep_ax_proc - mean(linep_ax_proc)).^2);
        Rsq = 1 - rn/sst;
        
        [~, se_ax] = nlparci2(a, real(res), 'Jacobian', real(J));
        
%         xC = a(3);
%         par_d = a(4);
%         par_e = a(5);
%         xC = a(2);
%         par_d = a(3);
%         par_e = a(4);
        xC = a(1);
        par_d = a(2);
        par_e = a(3);
        FWHM_ax = 2*sqrt(2)*par_d*(log(2))^(1/(2*par_e));
        x0ax = xC - FWHM_ax/2;
        x1 = xC + FWHM_ax/2;
        widthResult = [xC, x0ax, x1];
        % intensity = max(linep_ax);
        intensity = sum(linep_ax_proc);
        
        if param.plot_gauss
            tv = x(1):(x(2)-x(1))/10:x(end);
            plot(h_gauss, tv, superGaussian(a,tv))
        end
        
        %         [FWHM_ax, ~, ~, Rsq(ii), se_ax] = fitProfile_public(linep_ax,[],1, param); % axial FWHM
    end
    
%     fitvals_ax(ii,:) = [FWHM_ax se_ax(4)];
    fitvals_ax(ii,:) = [FWHM_ax se_ax(2)];
    
    %% Plot image with septal and axial axes
    
    if plot_im
        hold(h_im, 'on')
        
        slp = (y0(2)-y0(1)) / (x0(2)-x0(1));
        yint = y0(2) - (slp*x0(2));
        
        tv = min(x0):0.1:max(x0);
        plot(h_im, tv, tv*slp+yint, 'k', 'linew', 3)
        
        x_orth0 = [xy0(1)-r0*sin(theta) xy0(1)+r0*sin(theta)]; % for orthogonal axis straight through center of septum
        y_orth0 = [xy0(2)+r0*cos(theta) xy0(2)-r0*cos(theta)];
        
        slp_o = (y_orth0(2)-y_orth0(1)) / (x_orth0(2)-x_orth0(1));
        yint_o = y_orth0(2) - (slp_o*x_orth0(2));
        
        tv_o = min(x_orth0):0.1:max(x_orth0);
        plot(h_im, tv_o, tv_o*slp_o+yint_o, 'r', 'linew', 3)
        
        pause
    end
end
