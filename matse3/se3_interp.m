% interpolate 
% equivalent to slerp in pure quaternions
%
function y = se3_interp(x1,x2,alpha)

% exp(alpha*log(inv(x1)x2))*x1
y = se3_mul(se3_exp(alpha*se3_log(se3_mul(se3_inv(x1),x2))),x1)