% given a SE3 Gaussian evaluates the non linear function fx using Sigma Point.
% the output is another SE3
%
% g is the mean as se3
% S is the covariance 6x6
% f is the function
% sparams are the SP params (aka alpha,beta)
%
% This code can be modified to other cases such as:
% f: lie -> non-lie
% f: non-lie -> lie
%
% or generalized to multiple lie inputs
%
% Emanuele Ruffaldi
function [yd,xyS] = se3d_tx(d,f,sparams)

if nargin == 3
	sparams = [];
end

[g,S] = se3d_get(d);

[xp,wei] = se3d_sigmas(g,S,sparams);

% then evaluate the function
yp = zeros(size(xp,1),4,4);
for I=1:size(xp,1)
	yp(I,:,:) = f(squeeze(xp(I,:,:)));
end

% reassembly and compute co-variance
[yd,xyS] = se3d_unsigmas(yp,xp,wei);
