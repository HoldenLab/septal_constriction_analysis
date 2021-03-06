% Author: Kevin Whitley
% Date created: 200128

% This is an objective function to minimize for fitting
% horizontally-oriented bacterial septal profiles.

function obj = septum_1D_obj_cont(prof_im, guess, range, param)

if nargin == 0
    guess = [10 15 1 0];
    range = 1:length(prof_im);
end

mu = guess(1);
R = guess(2);
bg = guess(3);
amp = guess(4);
% ring_grad = guess(3);
psfFWHM = param.psfFWHM;
pixSz = param.pixSz;

prof_model = septum_model_1D_cont(mu, R, psfFWHM, pixSz, range, bg, amp);

sp = range(2)-range(1); % spacing
prof_model = downsample(prof_model, uint16(1/sp)); % same number of points as image profile

diff = prof_im - prof_model;

obj = sum(diff.^2);