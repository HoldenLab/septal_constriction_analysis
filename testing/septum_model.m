% Author: Kevin Whitley
% Date created: 190906

% This function provides a model of a 'rotated circle' representing a
% bacterial septum lying orthogonally to the imaging plane with center [mux
% muy], radius R, width W, and angle theta.

function im_model = septum_model(mu, R, W, theta, amp, ring_grad, psfFWHM, pixSz, range)

if nargin==0
    mu = [20 20];
    R = .5; % [pix]
    W = .5; % [pix] width of Z ring
    theta = pi/3; % note that -pi/4 is a discontinuity. seems to work fine anyway.
    amp = 1;
    ring_grad = 0;
    psfFWHM = 250;
    pixSz = 65;
    sp = 0.01;
    range = 1:sp:35;
end
% sp = range(2)-range(1);

X0 = mu(1);
Y0 = mu(2);

sigma = psfFWHM/2.35/pixSz;

[X, Y] = meshgrid(range, range);

% make axis 1 and axis 2 (rotated by theta, orthogonal) for both i and i+1.
% axis 1 is along septum, and axis 2 is along cell length. also need
% midpoints (i+1/2) for later.
ax1_i = (X(1:end-1,1:end-1)-X0)*cos(theta) + (Y(1:end-1,1:end-1)-Y0)*sin(theta); % axis along septum
ax1_ip1 = (X(2:end,2:end)-X0)*cos(theta) + (Y(2:end,2:end)-Y0)*sin(theta);
ax1_ipm = (ax1_i + ax1_ip1)./2; % between xi and x(i+1)
ax2_i = (X(1:end-1,1:end-1)-X0)*sin(theta) - (Y(1:end-1,1:end-1)-Y0)*cos(theta); % axis along cell length
ax2_ip1 = (X(2:end,2:end)-X0)*sin(theta) - (Y(2:end,2:end)-Y0)*cos(theta); % axis along cell length
ax2_ipm = (ax2_i + ax2_ip1)./2; % between yi and y(i+1)

% calculate arc lengths s(i) and s(i+1), then difference between them (ring
% intensity profile). asin will give complex numbers, but only want real
% part. anything outside real range will be cut later.
arcl_i = 2*R * real(asin(ax1_i/R)); % arc length out to axi
arcl_ip1 = 2*R * real(asin(ax1_ip1/R)); % arc length out to ax(i+1)
ring_prof = abs(arcl_ip1 - arcl_i); % Delta(arc length)

% cut image down based on radius and width of septum. this will remove any
% non-real parts of ring_prof as well.
ring_prof(abs(ax1_ipm)>R | abs(ax2_ipm)>W) = 0;

% change to real pixel size and blur based on psf FWHM.
ring_prof = imgaussfilt(ring_prof,sigma/sp);
ring_prof = imresize(ring_prof, [range(end) range(end)]);

ring_prof_scale = ring_prof/sum(ring_prof(:));
ring_prof_scale = amp*ring_prof_scale;

im_model = ring_prof_scale;

if nargin==0
    figure
    imagesc(im_model)
end