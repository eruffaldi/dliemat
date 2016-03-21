%
% Exponential from algebra to group
%
function y = se3_exp(x)

omega = x(1:3);
u = x(4:6);

theta = norm(omega);
if abs(theta) < 1e-10
    A = 1;
    B = 0.5;
    C = 1/6;
else
    A = sinc(theta);
    B = (1-cos(theta))/(theta^2);
    C = (1-A)/(theta^2);
end

S = skew(omega);
            
R = eye(3) + A*S + B*S^2;
V = eye(3) + B*S + C*S^2;

y = eye(4);
y(1:3,1:3) = R;
y(1:3,4) = V*u(:);            

function S = skew(v)
S = [  0   -v(3)  v(2)
      v(3)  0    -v(1)
     -v(2) v(1)   0];

