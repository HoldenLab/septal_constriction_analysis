% Author: Kevin Whitley
% Date created: 190724

% This objective function returns the square difference between an input
% image and a model ellipse top-hat.

function obj_var = bethfit_obj(im, guess, range)

if nargin == 0
    guess = [10 6 15 18 6 5e3];
    range = 1:31;
end

mu = [guess(1) guess(2); guess(3) guess(4)];
sigma = cat(3,[guess(5) guess(5)],[guess(5) guess(5)]);
ringbg = guess(6);

elrad1 = norm([mu(1,1) mu(1,2)] - [mu(2,1) mu(2,2)]) / 2; % major radius
elrad2 = sigma(1,1,1) /2; % minor radius
th = atan((mu(2,2)-mu(1,2)) / (mu(2,1)-mu(1,1)));
x0 = mu(1,1) + (mu(2,1)-mu(1,1))/2;
y0 = mu(1,2) + (mu(2,2)-mu(1,2))/2;

map = blurredEllipseTophat([x0 y0 elrad1 elrad2 th 1 ringbg], [length(range) length(range)]);

diff = im - map;

obj_var = sum(sum(diff.^2));