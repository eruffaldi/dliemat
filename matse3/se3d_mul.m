%
% Product
%
% Emanuele Ruffaldi SSSA 2015
function y=se3d_mul(a,b)

[ga,ca] = se3d_get(a);
[gb,cb] = se3d_get(b);

A = se3_adj(ga);
y = se3d_set(ga*gb,ca + A*cb*A'); % instead of: se3_mul(ga,gb)