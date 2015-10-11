% inverse of the distribution
function y=se3d_inv(a)

[ga,ca] = se3d_get(a); % extract 
A = se3_adj(gy);    % adjoint of inverse
y = se3d_set(se3_inv(ga),A*ca*A');
