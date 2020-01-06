% Author: Kevin Whitley
% Date created: 190722

% This function creates a map of a mixed Gaussian function with mean
% positions mu and covariance matrix sigma.

function pdfim = gmfun(mu, sigma, amp, ringbg, range)

if nargin==0
    mu = [10 6; 15 18];
    sigma = cat(3,[6 6],[6 6]);
    range = 1:31;
    amp = 5e3;
    ringbg = 10;
end

gm = gmdistribution(mu,sigma);

[X1, X2] = meshgrid(range, range);
X = [X1(:) X2(:)];
gmpdf = amp*pdf(gm, X);
gmim = reshape(gmpdf, [length(X1) length(X2)]);

elrad1 = norm([mu(1,1) mu(1,2)] - [mu(2,1) mu(2,2)]) / 2; % major radius (actually distance to focus)
elrad2 = sigma(1,1,1) /2; % minor radius
elrad_major = sqrt(elrad1^2 + elrad2^2); % true radius
th = atan((mu(2,2)-mu(1,2)) / (mu(2,1)-mu(1,1)));
x0 = mu(1,1) + (mu(2,1)-mu(1,1))/2;
y0 = mu(1,2) + (mu(2,2)-mu(1,2))/2;

ell = blurredEllipseTophat([x0 y0 elrad_major elrad2 th 0.5 ringbg], [length(range) length(range)]);

pdfim = gmim + ell;