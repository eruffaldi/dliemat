%
% Decomposes the distribution as mean and variance
%
function [g,S] = se3d_get(x)

g = reshape(x(1:16),4,4);
S = reshape(x(17:end),6,6);