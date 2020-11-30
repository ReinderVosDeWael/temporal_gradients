function y = eccentricity(x)
% Computes eccentricity of the first three vectors of x. 
%
%   ECCENTRICITY(x) computes the eccentricity i.e. square root of the sum
%   of squares of the first three vectors of x. x is an n-by-m gradient
%   matrix where n represents regions/vertices and m represents gradient
%   numbers.

y = sqrt(sum(x(:,1:3).^2,2));
end
