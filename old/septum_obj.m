% Author: Kevin Whitley
% Date created: 190909

% This is an objective function to minimize for fitting
% horizontally-oriented bacterial septa.

function obj = septum_obj(im, guess, range, param)

if nargin == 0
    guess = [10 6 15 18 6 1 0.2 5e3 1];
    range = 1:31;
end

mu = [guess(1) guess(2)];
R = guess(3);
W = guess(4);
theta = guess(5);
amp = guess(6);
ring_grad = guess(7);
psfFWHM = param.psfFWHM;
pixSz = param.pixSz;

map = septum_model(mu, R, W, theta, amp, ring_grad, psfFWHM, pixSz, range);

diff = im - map;

obj = sum(sum(diff.^2));