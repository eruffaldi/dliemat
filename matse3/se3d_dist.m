% Estimate distance between two SE3 gaussians
%
%
% Emanuele Ruffaldi 2016
function [d,S] = se3d_dist(a,b,order)

if nargin == 2
    order = 2;
end
ab = se3d_mul(sed3d_inv(a),b,order);

mu_ab = se3d_mean(ab);
S = se3d_cov(ab);
d = se3_log(mu_ab);
