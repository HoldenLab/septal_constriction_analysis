% Author: Kevin Whitley
% Date created: 190722

% This objective function returns the square difference between an input
% image and a model mixed Gaussian function (2 Gaussians only).

function obj_var = gmfit_obj(im, guess, range)

if nargin == 0
    guess = [10 6 15 18 6 1 5e3];
    range = 1:31;
end

mu = [guess(1) guess(2); guess(3) guess(4)];
sigma = cat(3,[guess(5) guess(5)],[guess(5) guess(5)]);
amp = exp(guess(6));
ringbg = guess(7);

map = gmfun(mu, sigma, amp, ringbg, range);

diff = im - map;

obj_var = sum(sum(diff.^2));