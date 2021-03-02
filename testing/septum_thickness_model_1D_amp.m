% Author: Kevin Whitley
% Date created: 200128

% This function produces a line profile of a septum using an explicit
% 'tilted circle' model.

function improf = septum_thickness_model_1D_amp(X0, R, psfFWHM, pixSz, X, bg, amp)

if nargin==0
    X0 = 19.5335;
    R = 3.1; % [pix]
    psfFWHM = 250;
    pixSz = 65;
    bg = 2.3;
    amp = 12.4;
    sp = .01; % spacing needs to be low enough to prevent problems with discretization (probably <0.5 or so)
    X = 1:sp:35;
end

sp = X(2)-X(1);

sigma = psfFWHM/2.35/pixSz;

thick_prof = zeros(1, length(X)-1);
thickline = ones(1, round(2*R/sp));

pt1 = max([1 round(X0/sp-length(thickline)/2+1)]);
pt2 = min([round(X0/sp+length(thickline)/2) length(X)-1]);
thick_prof(pt1:pt2) = ones(1, length(pt1:pt2));

gpdf_fcn = @(a,x) 1/(a(2)*sqrt(2*pi)) * exp(-(x-a(1)).^2 ./ (2*a(2)^2));
gauss = gpdf_fcn([0 sigma], X-X0);

thick_prof = conv(thick_prof, gauss, 'same'); % convolve with Gaussian psf

ring_prof_scale = thick_prof/max(thick_prof(:));

improf = ring_prof_scale * amp + bg;

if nargin==0
    figure
    hold on
    plot(X(1:end-1), improf)
end