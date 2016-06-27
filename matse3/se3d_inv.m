% inverse of the distribution
%
% Emanuele Ruffaldi 2015
function y=se3d_inv(a)

[ga,S] = se3d_get(a); % extract 
iga = se3_inv(ga);
A = se3_adj(iga);    % adjoint of inverse
y = se3d_set(iga,A*S*A');
