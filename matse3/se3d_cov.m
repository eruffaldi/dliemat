function [S] = se3d_cov(x)

S = reshape(x(17:end),6,6);