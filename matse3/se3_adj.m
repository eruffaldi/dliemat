%
% Computes the Adjoint
%
function y = se3_adj(x)

R = x(1:3,1:3);
t = x(1:3,4);
% SKEW for [t,omega] y = [R skew(t)*R; zeros(3) R];
% SKEW for [omega,t] y = [R zeros(3);skew(t)*R R];

y = [R skew(t)*R; zeros(3) R];


function S = skew(v)
S = [  0   -v(3)  v(2)
      v(3)  0    -v(1)
     -v(2) v(1)   0];
