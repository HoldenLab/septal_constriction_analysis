% Author: Kevin Whitley
% Date created: 200128

% This function produces a line profile of a septum using an explicit
% 'tilted circle' model.

function improf = septum_model_1D(X0, R, psfFWHM, pixSz, X)

if nargin==0
    X0 = 19.5335;
    R = 8.7665; % [pix]
    ring_grad = 0;
    psfFWHM = 250;
    pixSz = 65;
    sp = .01; % spacing needs to be low enough to prevent problems with discretization (probably <0.5 or so)
    X = 1:sp:25;
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

improf = ring_prof_scale;

if nargin==0
%     figure
    hold on
    halfX = X - (X(2)-X(1))/2;
    plot(halfX, improf)
end