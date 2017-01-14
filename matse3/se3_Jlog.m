% 3 x 12
function [v,J] = se3_Jlog(A)

v = se3_log(A);
J = vec2jac(v); % 3x3