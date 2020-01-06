% Author: Kevin Whitley
% Date created: 190723

% This function fits an image of a bacterial septum to an explicit model
% using a minimization routine.

function [fin, cell_or, fval] = fit_septum_explicit(imstack, plot_im, param)

fval_thresh = param.fval_thresh;
ssim_thresh = param.ssim_thresh;
cc_thresh = param.cc_thresh;

psfFWHM = param.psfFWHM; % [nm]
pixSz = param.pixSz; % [nm/pixel]
gridsp = param.gridsp;

gw = psfFWHM/2.35/pixSz; % Gaussian width (sigma)

if plot_im
    figure
    h_im = gca;
end

fin = []; cell_or = []; sept_or = []; ssimval = [];
for ii = 1:size(imstack,3)
    
    frame = imstack(:,:,ii);
    
    if plot_im
        hold off
        imagesc(h_im, frame)
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
        if cc2.NumObjects == 2 % only change to new bwframe if there are now two objects
            bwframe = bwframe2;
        end
    end
    
    % remove all but one or two blobs from image
    BW = bwareafilt(bwframe, 2, 'largest');
    fbgsub(BW==0) = 0;
    bwframe(BW==0) = 0;
    
    % are there really two objects? if one is way bigger than the other
    % it's probably not a good septum image. should filter out some crap.
    cc = bwconncomp(bwframe);
    if cc.NumObjects == 2
        props = regionprops(bwframe, 'Area');
        area_rat = max(props.Area)/min(props.Area);
        if area_rat > 10
            BW = bwareafilt(bwframe, 1, 'largest');
            fbgsub(BW==0) = 0;
            bwframe(BW==0) = 0;
        end
    end
    
    % Deprecated 190913
%     bwframe2 = seg_f > 1;
%     bwframe2_or = regionprops(bwframe2,'Orientation');
%     cell_or(ii) = max(cat(1,bwframe2_or.Orientation)); % get orientation of cell by cytoplasmic fluorescence
    
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
            continue
        end
        
        if param.gm_model            
            b0 = [cents(1,:) cents(2,:) gw 0.5 0.5 0];
         else
            cent_dist = (cents(2,:) - cents(1,:)) / 2; % distance from centroid to centroid
            cent_pos = cents(2,:) - cent_dist; % position of center between centroids
            guess_th = atan((cents(2,2)-cents(1,2))./(cents(2,1)-cents(1,1))); % orientation of septum
            guess_r = norm(cents(1,:)-cents(2,:))/2;
        end
    elseif size(cents,1) == 1 % mature ring
        
        % Deprecated 190913
        %         sept_or(ii) = max(cat(1, props.Orientation));

        if param.gm_model
            [yf, xf] = find(bwframe);
            edge1 = [xf(1) yf(1)];
            edge2 = [xf(end) yf(end)];
            guessp1 = (edge1+cents)/2;
            guessp2 = (edge2+cents)/2;
            
            b0 = [guessp1 guessp2 gw 0.5 0.3 0];
        else
            cent_pos = cents;
            guess_th = -props.Orientation * pi/180; % minus sign because it needs to be reflected about y axis.
            %             guess_th = atan((edge2(2)-edge1(2))/(edge2(1)-edge1(1)));
            guess_r = props.MajorAxisLength/2 * 0.7;
            %             guess_r = norm(edge2-edge1)/2*0.7;
            if isnan(guess_th) && ii~=1
                guess_th = fin(ii-1,4);
                guess_r = 0.2;
            end
        end
    else % probably no centroid because image is utter crap
        fin(ii,:) = ones(1,4)*NaN;
        continue
    end
    
    %% Fit
    
    s_frame = size(fbgsub,1);
    im_scaled = fbgsub/sum(fbgsub(:));

    if param.gm_model
        lb = [1 1 1 1 gw 0 0.2 0];
        ub = [31 31 31 31 gw 5 0.8 1];
        options = optimoptions('fmincon', 'OptimalityTolerance', 1e-11, 'MaxFunctionEvaluations', 1e3, 'StepTolerance', 1e-20, 'Display', 'Iter-Detailed');
        
        [fin(ii,:), fval(ii)] = fmincon(@(x)gmfit_obj(im_scaled,x,1:gridsp:s_frame), b0, [], [], [], [], lb, ub, [], options);
        mod = gmfun([fin(ii,1) fin(ii,2); fin(ii,3) fin(ii,4)], cat(3,[fin(ii,5) fin(ii,5)],[fin(ii,5) fin(ii,5)]), exp(fin(ii,6)), fin(ii,7), fin(ii,8), 1:gridsp:s_frame);
        gues = gmfun([b0(1) b0(2); b0(3) b0(4)], cat(3,[b0(5) b0(5)],[b0(5) b0(5)]), exp(b0(6)), b0(7), b0(8), 1:gridsp:s_frame);
        
        if plot_im
            hold on
            plot(h_im, fin(ii,1), fin(ii,2), 'om')
            plot(h_im, fin(ii,3), fin(ii,4), 'ok')
        end
    else
        wid = 1.2;
        b0 = [cent_pos guess_r wid guess_th 1 0];
%         lb = [1 1 0.1 1.5 0 0.5 -1];
%         ub = [31 31 20 3 pi 10 1];
        
%         [fin(ii,:), fval(ii), ~, out(ii)] = fmincon(@(x)septum_obj(im_scaled,x,1:0.1:s_frame), b0, [], [], [], [], lb, ub, [], options);
        [fin(ii,:), fval(ii)] = fminsearch(@(x)septum_obj(im_scaled,[x(1:3) wid x(4) 1 0],1:gridsp:s_frame, param), [b0(1:3) b0(5)]);
        mod = septum_model([fin(ii,1) fin(ii,2)], fin(ii,3), wid, fin(ii,4), 1, 0, psfFWHM, pixSz, 1:gridsp:s_frame);
        gues = septum_model([b0(1) b0(2)], b0(3), b0(4), b0(5), b0(6), b0(7), psfFWHM, pixSz, 1:gridsp:s_frame); % for debugging
        
        if plot_im
            hold on
            plot(h_im, fin(ii,1)+fin(ii,3)*cos(fin(ii,4)), fin(ii,2)+fin(ii,3)*sin(fin(ii,4)), 'xk', 'linew', 2, 'MarkerSize', 12)
            plot(h_im, fin(ii,1)-fin(ii,3)*cos(fin(ii,4)), fin(ii,2)-fin(ii,3)*sin(fin(ii,4)), 'xk', 'linew', 2, 'MarkerSize', 12)
        end
    end

    %% Filter out crap
    
    % filter out poor fits
    if fval(ii) > fval_thresh
        fin(ii,:) = ones(1,size(fin(ii,:),1))*NaN;
    end
    
    % filter based on structural similarity matrix (compare to bg
    % subtracted image)
    modnorm = mod/sum(mod(:));
    fnorm = fbgsub2/sum(fbgsub2(:));
    [ssimval(ii), ~] = ssim(fnorm, modnorm);
    if ssimval(ii) < ssim_thresh
        fin(ii,:) = ones(1,size(fin(ii,:),1))*NaN;
    end
    
    % filter based on cross-correlation between original image and model
    normccmap = normxcorr2(frame/sum(frame(:)), modnorm);
    hlfnrm = floor(size(normccmap,2)/2);
    nrmcutsz = 4; % [pix]
    cutnormmap = normccmap(hlfnrm-nrmcutsz:hlfnrm+nrmcutsz, hlfnrm-nrmcutsz:hlfnrm+nrmcutsz);
    ccval(ii) = max(cutnormmap(:));
    if ccval(ii) < cc_thresh
        fin(ii,:) = ones(1,size(fin(ii,:),1))*NaN;
    end
    
end
