% Author: Kevin Whitley
% Date: 190919

% This function fits constriction vs. time data to a choice model.

function [fits, rn, se, fcn] = fit_constriction_data(tdat, ydat, model)

switch model
    case 'linear'
        % fcn = @(a,t) a(3).*heaviside(a(1)-t) + heaviside(t-a(1)).*(a(3)+a(2)*(t-a(1)));
        fcn = @(a,t) ydat(1).*heaviside(a(1)-t) + heaviside(t-a(1)).*(ydat(1)+a(2)*(t-a(1))); % 2 param fit
        b0 = [tdat(1) -40];
        lb = [0 -100];
        ub = [tdat(end) 0];
    case 'parabolic'
        fcn = @(a,t) a(3).*heaviside(a(1)-t) + real(heaviside(t-a(1)).*sqrt(a(3)^2 - a(2)*(t-a(1))));
        b0 = [(tdat(end)-tdat(1))/2 5e4 900];
%         lb = [tdat(1) 5e3 900];
        lb = [tdat(1)-5 5e3 900];
        ub = [tdat(end) 3e5 1300];
    case 'logistic'
        nt = @(a,t) (t + 1/a(1)*log(abs((1+exp(-a(1)*(t-a(2))))/(1+exp(a(1)*a(2)))))); % total synthase activity over time (integral of generic logistic function)
        %                 fcn = @(a,t) sqrt(w0^2 - a(1)*nt([a(2) a(3)],t));
        fcn = @(a,t) real(sqrt(a(4)^2 - a(2)*nt([a(3) a(1)],t)));
        b0 = [(tdat(end)-tdat(1))/2 5e4 1 ydat(1)];
        lb = [tdat(1) 5e3 1 700];
        ub = [tdat(end) 3e5 10 1200];
    case 'kludge_model'
        fcn = @(a,t) a(3).*heaviside(a(1)-t) + real(heaviside(t-a(1)).*(a(3)^(2+a(4)) - a(2)*(t-a(1))).^(1/(2+a(4))));
        b0 = [(tdat(end)-tdat(1))/2 5e4 ydat(1) 1];
        lb = [tdat(1) 5e3 700 0];
        ub = [tdat(end) 3e5 1200 10];
end

options = optimoptions('lsqcurvefit','Display','off');   
[fits, rn, res, ~, ~, ~, J] = lsqcurvefit(fcn, b0, tdat, ydat, lb, ub, options);
[~, se] = nlparci2(fits, real(res), 'Jacobian', real(J));