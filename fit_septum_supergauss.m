% Author: Kevin Whitley
% Date created: 190723

% This function fits an image of a bacterial septum to an explicit model
% using a minimization routine.

function [improf, orthprof, fitex] = fit_septum_supergauss(imstack, plot_im, param)

% fval_thresh = param.fval_thresh;
% ssim_thresh = param.ssim_thresh;
% cc_thresh = param.cc_thresh;

% psfFWHM = param.psfFWHM; % [nm]
% pixSz = param.pixSz; % [nm/pixel]
% gridsp = param.gridsp;

% gw = psfFWHM/2.35/pixSz; % Gaussian width (sigma)

% plot_gauss = param.plot_gauss;
plot_gauss = 0;

if plot_im
    figure
    h_im = gca;
end
if plot_gauss
    figure
    h_gauss = gca;
end

fin = []; ssimval = []; improf = []; orthprof = [];
for ii = 1:size(imstack,3)
    
    frame = imstack(:,:,ii);
    
    if plot_im
        hold(h_im, 'off')
        imagesc(h_im, frame)
    end    
    if plot_gauss
        hold(h_gauss, 'off')
    end
    
    %% Process image before fitting - try to find two blobs
    
    % custom thresh
    %     thr = max(max(frame)) / 1.7; % threshold
    %     bwframe = frame > thr; % binary image
    %     fbgsub = frame - thr;
    %     fbgsub(~bwframe) = 0;
    
    % otsu thresh
    %     T = otsu(frame);
    %     fbgsub = frame;
    %     fbgsub(T==1) = 0;
    %     bwframe = fbgsub > 0;
    
    % multilevel thresh
    thr = multithresh(frame, 3);
    seg_f = imquantize(frame, thr);
    
    % background subtracted image
    fbgsub = frame;
    fbgsub(seg_f<=2) = 0;
    
    % less background subtracted image (used for assessing image quality
    % after fitting)
    fbgsub2 = frame;
    fbgsub2(seg_f<=1) = 0;
    
    bwframe = seg_f > length(thr); % top level
    
    % if bwframe only has one object, drop down one level and try again (helps with asymmetric nascent septa)
    cc = bwconncomp(bwframe);
    if cc.NumObjects < 2
        bwframe2 = seg_f > length(thr)-1;
        cc2 = bwconncomp(bwframe2);
%         if cc2.NumObjects == 2 % only change to new bwframe if there are now two objects
            bwframe = bwframe2;
%         end
    end
    
    % remove all but one or two blobs from image
    BW = bwareafilt(bwframe, 2, 'largest');
    fbgsub(BW==0) = 0;
    bwframe(BW==0) = 0;
    
    % are there really two objects? if one is way bigger than the other
    % it's probably not a good septum image. should filter out some crap.
    
    % Deprecated 191211
%     cc = bwconncomp(bwframe);
%     if cc.NumObjects == 2
%         props = regionprops(bwframe, 'Area');
%         area_rat = max(props.Area)/min(props.Area);
%         if area_rat > 10
%             BW = bwareafilt(bwframe, 1, 'largest');
%             fbgsub(BW==0) = 0;
%             bwframe(BW==0) = 0;
%         end
%     end
    
    % blur
    fbgsub = imgaussfilt(fbgsub,1.5);
    fbgsub2 = imgaussfilt(fbgsub2, 1.5);
    
    %% Calculate guess values before fitting

    props = regionprops(bwframe, 'Centroid', 'Area', 'Orientation', 'MinorAxisLength', 'MajorAxisLength', 'Area');
    cents = cat(1, props.Centroid); % centroids in image
    
    % the fits are quite sensitive to initial guesses, so choose these
    % carefully
    if size(cents,1) == 2 % nascent ring
        
        if props(1).Area>65 || props(2).Area>65
            fin(ii,:) = ones(1,4)*NaN;
            fval(ii) = NaN;
            improf(ii,:) = ones(1,size(improf,2))*NaN;
            continue
        end
        cent_dist = (cents(2,:) - cents(1,:)) / 2; % distance from centroid to centroid
        cent_pos = cents(2,:) - cent_dist; % position of center between centroids
        guess_th = atan((cents(2,2)-cents(1,2))./(cents(2,1)-cents(1,1))); % orientation of septum
        guess_r = norm(cents(1,:)-cents(2,:))/2;
        
    elseif size(cents,1) == 1 % mature ring
        
        cent_pos = cents;
        guess_th = -props.Orientation * pi/180; % minus sign because it needs to be reflected about y axis.
        %             guess_th = atan((edge2(2)-edge1(2))/(edge2(1)-edge1(1)));
        guess_r = props.MajorAxisLength/2 * 0.7;
        %             guess_r = norm(edge2-edge1)/2*0.7;
        if isnan(guess_th) && ii~=1
            guess_th = fin(ii-1,4);
            guess_r = 0.2;
        end
    else % probably no centroid because image is utter crap
        fin(ii,:) = ones(1,4)*NaN;
        improf(ii,:) = ones(1,size(improf,2))*NaN;
        continue
    end
    
    %% Fit
    
    s_frame = size(fbgsub,1);
    im_scaled = fbgsub/sum(fbgsub(:));
    
    wid = 1.2;
    b0 = [cent_pos guess_r wid guess_th 1 0];
    
    r0 = 1800 / 2 / 65; % [pix]
    x0(ii,:) = [b0(1)+r0*cos(b0(5)) b0(1)-r0*cos(b0(5))];
    y0(ii,:) = [b0(2)+r0*sin(b0(5)) b0(2)-r0*sin(b0(5))];
    
    x_orth(ii,:) = [b0(1)-r0*sin(b0(5)) b0(1)+r0*sin(b0(5))];
    y_orth(ii,:) = [b0(2)+r0*cos(b0(5)) b0(2)-r0*cos(b0(5))];
    %         lb = [1 1 0.1 1.5 0 0.5 -1];
    %         ub = [31 31 20 3 pi 10 1];
    
    %         [fin(ii,:), fval(ii), ~, out(ii)] = fmincon(@(x)septum_obj(im_scaled,x,1:0.1:s_frame), b0, [], [], [], [], lb, ub, [], options);
%     [fin(ii,:), fval(ii)] = fminsearch(@(x)septum_obj(im_scaled,[x(1:3) wid x(4) 1 0],1:gridsp:s_frame, param), [b0(1:3) b0(5)]);
%     mod = septum_model([fin(ii,1) fin(ii,2)], fin(ii,3), wid, fin(ii,4), 1, 0, psfFWHM, pixSz, 1:gridsp:s_frame);
%     gues = septum_model([b0(1) b0(2)], b0(3), b0(4), b0(5), b0(6), b0(7), psfFWHM, pixSz, 1:gridsp:s_frame); % for debugging

%     x0(ii,:) = [fin(ii,1)+r0*cos(fin(ii,4)) fin(ii,1)-r0*cos(fin(ii,4))];
%     y0(ii,:) = [fin(ii,2)+r0*sin(fin(ii,4)) fin(ii,2)-r0*sin(fin(ii,4))];
    
    ip = improfile(frame, x0(ii,:), y0(ii,:))';
    
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
    
    imp = improf(ii,~isnan(improf(ii,:)));
    imp_proc = (imp-min(imp));
    imp_proc = imp_proc/max(imp_proc);
    [pks, locs] = findpeaks(imp_proc);
    pkmat = [pks' locs'];
    pkmat = sortrows(pkmat);
    if size(pkmat,1)>1
        initwidth = abs((pkmat(end)-pkmat(end-1))/1.5);
    else
        initwidth = 3;
    end
    initguess = [length(imp)/2+1 initwidth];
    options = optimset('Display', 'off');
    range = 1:length(imp);
%     lb = [0 0];
%     ub = [length(imp) 15];
%     range = -length(imp)/20:0.1:length(imp)/20 - 0.1;
%     fitex = lsqcurvefit(@(a,x)septum_model_1D(a(1),a(2),250,65,x), initguess, range, imp_proc, lb, ub);
    fitex(ii,:) = fminsearch(@(x)septum_1D_obj(imp_proc,x,range,param), initguess, options);
%     fitex = fmincon(@(x)septum_1D_obj(imp_proc,x,range,param), initguess, [], [], [], [], lb, ub);
    
    if param.plot_explicit
        figure
        hold on
        plot(range, imp_proc)
%         tr = 1:0.1:length(imp);
%         plot(tr,septum_model_1D(fitex(1),fitex(2),250,65,tr))
        plot(range, septum_model_1D(fitex(ii,1),fitex(ii,2),250,65,range))
    end
    
    % Get line profiles along cell axis to get width of septum in that
    % direction
    slp = (y0(ii,2)-y0(ii,1)) / (x0(ii,2)-x0(ii,1));
    yint = y0(ii,2) - (slp*x0(ii,2));
    if isinf(slp)
        orthprof(ii,:) = ones(size(orthprof,2),1)*NaN;
        continue
    end
    orth_space = 0.1;
    x0_orth_ind = min([x0(ii,1) x0(ii,2)]):orth_space:max([x0(ii,1) x0(ii,2)]);
    x_orth2 = []; y_orth2 = []; ip_orth = [];
    for jj = 1:length(x0_orth_ind)
        
        % new 'central' coordinates (along septum axis)
        new_x0 = x0_orth_ind(jj);
        new_y0 = new_x0*slp + yint;
        
        x_orth2(jj,:) = [new_x0-r0*sin(b0(5)) new_x0+r0*sin(b0(5))];
        y_orth2(jj,:) = [new_y0+r0*cos(b0(5)) new_y0-r0*cos(b0(5))];

        % line profile along cell axis with shifted 'central' coordinates
        ip_orth(jj,:) = improfile(frame, x_orth2(jj,:), y_orth2(jj,:))';

        if plot_im && 0 % for debugging
            hold(h_im, 'on')
            slp_o1 = (y_orth2(jj,2)-y_orth2(jj,1)) / (x_orth2(jj,2)-x_orth2(jj,1));
            yint_o1 = y_orth2(jj,2) - (slp_o1*x_orth2(jj,2));
            tv_o1 = min([x_orth2(jj,1) x_orth2(jj,2)]):max([x_orth2(jj,1) x_orth2(jj,2)]);
            plot(h_im, tv_o1, tv_o1*slp_o1+yint_o1, 'r', 'linew', 3)
        end
        if plot_gauss && 0
            hold(h_gauss, 'on')
            plot(h_gauss, ip_orth(jj,:))
        end
    end
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
    orthprof(ii,:) = nansum(ip_orth,1)*orth_space; % integral rather than sum. will not give accurate intensity measurement!

    if plot_im
        hold(h_im, 'on')
%         plot(h_im, fin(ii,1)+fin(ii,3)*cos(fin(ii,4)), fin(ii,2)+fin(ii,3)*sin(fin(ii,4)), 'xk', 'linew', 2, 'MarkerSize', 12)
%         plot(h_im, fin(ii,1)-fin(ii,3)*cos(fin(ii,4)), fin(ii,2)-fin(ii,3)*sin(fin(ii,4)), 'xk', 'linew', 2, 'MarkerSize', 12)
        
%         plot(h_im, x0(ii,1),y0(ii,1),'xk', 'linew', 2, 'MarkerSize', 12)
%         plot(h_im, x0(ii,2),y0(ii,2),'xk', 'linew', 2, 'MarkerSize', 12)
        
%         plot(h_im, x_orth(ii,1),y_orth(ii,1),'ok', 'linew', 2, 'MarkerSize', 12)
%         plot(h_im, x_orth(ii,2),y_orth(ii,2),'ok', 'linew', 2, 'MarkerSize', 12)
        
        slp = (y0(ii,2)-y0(ii,1)) / (x0(ii,2)-x0(ii,1));
        yint = y0(ii,2) - (slp*x0(ii,2));
        
        tv = min([x0(ii,1) x0(ii,2)]):max([x0(ii,1) x0(ii,2)]);
        plot(h_im, tv, tv*slp+yint, 'k', 'linew', 3)
        
        slp_o = (y_orth(ii,2)-y_orth(ii,1)) / (x_orth(ii,2)-x_orth(ii,1));
        yint_o = y_orth(ii,2) - (slp_o*x_orth(ii,2));
        
        tv_o = min([x_orth(ii,1) x_orth(ii,2)]):max([x_orth(ii,1) x_orth(ii,2)]);
        plot(h_im, tv_o, tv_o*slp_o+yint_o, 'r', 'linew', 3)
    end
    if plot_gauss
        plot(h_gauss, improf(ii,:))
    end
    
    %% Filter out crap
    
    % filter out poor fits
%     if fval(ii) > fval_thresh
%         fin(ii,:) = ones(1,size(fin(ii,:),1))*NaN;
%     end
%     
%     % filter based on structural similarity matrix (compare to bg
%     % subtracted image)
%     modnorm = mod/sum(mod(:));
%     fnorm = fbgsub2/sum(fbgsub2(:));
%     [ssimval(ii), ~] = ssim(fnorm, modnorm);
%     if ssimval(ii) < ssim_thresh
%         fin(ii,:) = ones(1,size(fin(ii,:),1))*NaN;
%     end
%     
%     % filter based on cross-correlation between original image and model
%     normccmap = normxcorr2(frame/sum(frame(:)), modnorm);
%     hlfnrm = floor(size(normccmap,2)/2);
%     nrmcutsz = 4; % [pix]
%     cutnormmap = normccmap(hlfnrm-nrmcutsz:hlfnrm+nrmcutsz, hlfnrm-nrmcutsz:hlfnrm+nrmcutsz);
%     ccval(ii) = max(cutnormmap(:));
%     if ccval(ii) < cc_thresh
%         fin(ii,:) = ones(1,size(fin(ii,:),1))*NaN;
%     end
    
end
