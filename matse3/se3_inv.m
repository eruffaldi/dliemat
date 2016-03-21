% Inverse of the group
%
function y = se3_inv(x)

R = x(1:3,1:3);
t = x(1:3,4);


y = eye(4);
y(1:3,1:3) = R';
y(1:3,4) = -R'*t;