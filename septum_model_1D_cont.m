% Author: Kevin Whitley
% Date created: 200128

% This function produces a line profile of a septum using an explicit
% 'tilted circle' model.

function improf = septum_model_1D_cont(X0, R, psfFWHM, pixSz, X, bg, amp)

if nargin==0
    X0 = 20.0049;
    R = 6.2651; % [pix]
    ring_grad = 0;
    psfFWHM = 250;
    pixSz = 65;
    bg = 0;
    amp = 1;
    sp = 0.001; % spacing needs to be low enough to prevent problems with discretization (probably <0.5 or so)
    X = 1:sp:39;
end

sigma = psfFWHM/2.35/pixSz;

ring_prof = 1./sqrt(1-((X-X0)/R).^2);

ring_prof = real(ring_prof);
ring_prof(isinf(ring_prof)) = 0;

gpdf_fcn = @(a,x) 1/(a(2)*sqrt(2*pi)) * exp(-(x-a(1)).^2 ./ (2*a(2)^2));
gauss = gpdf_fcn([0 sigma], X-X0);

ring_prof = conv(ring_prof, gauss, 'same'); % convolve with Gaussian psf
% ring_prof = ring_prof(floor(length(gauss)/2):floor((end-length(gauss)/2))+1); % cut ends (convolution increases length of vector)

ring_prof_scale = ring_prof/max(ring_prof(:));

improf = ring_prof_scale * amp + bg;

if nargin==0
    figure
    hold on
    plot(X, improf)
end