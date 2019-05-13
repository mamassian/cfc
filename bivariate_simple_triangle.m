function VV = bivariate_simple_triangle(hh, kk)

% Compute the integral of the independent (rho = 0) bivariate Normal
%   distribution with parameters (mx = 0, my = 0, sigma = 1) over
%   the right triangle bounded by (y = 0), (x = h), and (y = x*k/h).
%
% Based on the formula in Abramowitz & Stegun (1972), page 940.
%
% Uses the 'bvnu' function from Alan Genz (1998).
% Alternatively, use the 'normcdf' and 'mvncdf' functions from the
%   Statistics toolbox.
% 
% 08-AUG-2017 - created - pascal mamassian

aa = kk / hh;
rr = - aa / sqrt(1 + aa*aa);

VV = 1/4 + LL(hh, 0, rr) - LL(0, 0, rr) - (1/2)*QQ(hh);


% -> nested function for univariate cumulative up to +infinity
    function QQ_val = QQ(hh)
%         QQ_val = normcdf(-hh, 0, 1);
        QQ_val = (1 - erf(hh / sqrt(2))) / 2;
    end

% -> nested function for bivariate cumulative up to +infinity on x and y
    function LL_val = LL(hh, kk, cc)
%         LL_val = mvncdf([-hh, -kk], [0, 0], [1, cc; cc, 1]);
        LL_val = bvnu(hh, kk, cc);
    end

end
% *** THE END ***