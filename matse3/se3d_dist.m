% Estimate distance between two SE3 gaussians
%
%
% Emanuele Ruffaldi 2016
function [d,S] = se3d_dist(a,b)

ab = se3d_mul(sed3d_inv(a),b)

mu_ab = se3d_mean(ab)
S = se3d_cov(ab)
d = se3_log(mu_ab)
