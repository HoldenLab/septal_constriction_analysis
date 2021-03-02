function [FWHM, intensity, widthResult, Rsq, se] = fitProfile_public(im, width, plotpos_y, param)

if nargin<4
    param.plot_gauss = 0;
end

y = double(im(plotpos_y, :));
h = size(im, 2);
x = 1:h;

if param.plot_gauss
    figure
    h_gauss = gca;
    hold(h_gauss, 'on')
    plot(y)
end

superGaussian = @(a,x)  a(1) + a(2)*exp(-((x-a(3)).^2 / (2*a(4)^2)).^ a(5));
initwidth = 5;
initguess = [min(y), max(y), mean(x), initwidth, 1];
options = optimoptions('lsqcurvefit', 'Display', 'off');
lb = [0 0 0 0 0.5];
[a, rn, res, ~, ~, ~, J] = lsqcurvefit(superGaussian, initguess, x, y, lb, [], options);

sst = sum((y - mean(y)).^2);
Rsq = 1 - rn/sst;

[~, se] = nlparci2(a, real(res), 'Jacobian', real(J));

xC = a(3);
par_d = a(4);
par_e = a(5);
FWHM = 2*sqrt(2)*par_d*(log(2))^(1/(2*par_e));
x0 = xC - FWHM/2;
x1 = xC + FWHM/2;
widthResult = [xC, x0, x1];
% intensity = max(y);
intensity = sum(y);

if param.plot_gauss
    tv = x(1):(x(2)-x(1))/10:x(end);
    plot(h_gauss, tv, superGaussian(a,tv))
end

end

