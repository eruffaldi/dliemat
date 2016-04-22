%
% Exponential from algebra to group
%
function y = se3_exp(x)

omega = x(1:3);
u = x(4:6);

% first log of R
%simplify(taylor(sin(x)/x,x,0,'Order',6))
%       x^4/120 - x^2/6 + 1
%simplify(taylor((1-cos(x))/x/x,x,0,'Order',6))
%       x^4/720 - x^2/24 + 1/2
%simplify(taylor((1-sin(x)/x)/x/x,x,0,'Order',6))
%       x^4/5040 - x^2/120 + 1/6
theta = norm(omega);
if abs(theta) < 1e-10
    A = 1;
    B = 0.5;
    C = 1/6;
    S = zeros(3);
    R = eye(3) + A*S + B*S^2;
else
    A = sin(theta)/theta;
    B = (1-cos(theta))/(theta^2);
    C = (theta-sin(theta))/(theta^3);
    S = skew(omega);
    R = eye(3) + A*S + B*S^2;
    % TROUBLE IF NOT UNITARY
end

            
V = eye(3) + B*S + C*S^2;

y = eye(4);
y(1:3,1:3) = R;
y(1:3,4) = V*u(:);            

function S = skew(v)
S = [  0   -v(3)  v(2)
      v(3)  0    -v(1)
     -v(2) v(1)   0];

