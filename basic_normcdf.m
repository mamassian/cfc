% basic normcdf function with no checks, should be quicker than the
% built-in version.
% written by TB 23.05.20

function p = basic_normcdf(x,mu,sigma);

if nargin < 2
    mu = 0;
end
if nargin < 3
    sigma = 1;
end

z = (x-mu) ./ sigma;
p = 0.5 * erfc(-z ./ sqrt(2));

end