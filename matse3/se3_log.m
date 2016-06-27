%
% Logarithm of the group: the algebra
%
% Emanuele Ruffaldi 2016
%
function y = se3_log(x)

R = x(1:3,1:3);
t = x(1:3,4);

% first log of R
%simplify(taylor(sin(x)/x,x,0,'Order',6))
%       x^4/120 - x^2/6 + 1
%simplify(taylor((1-cos(x))/x/x,x,0,'Order',6))
%       x^4/720 - x^2/24 + 1/2
%simplify(taylor((1-sin(x)/x)/x/x,x,0,'Order',6))
%       x^4/5040 - x^2/120 + 1/6
theta = acos(max(-1,min((trace(R)-1)/2,1)));
assert(isreal(theta),'thetareal');

if abs(theta) < 1e-10
    B = 0.5;
    SO = (1/(2))*(R-R');  % =skew(omega)
    iV = eye(3); % use directly skew of omega
else
    A = sin(theta)/theta;
    B = (1-cos(theta))/(theta*theta);
    SO = (1/(2*A))*(R-R');  % =skew(omega)
    %??
    % syms x real
    % A = sin(x)/x
    % B= (1-cos(x))/(x*x)
    % Q=1/(x*x)*(1 - A/2/B)
    % taylor(Q,x,0)
    %       x^4/30240 + x^2/720 + 1/12
    Q= 1/(theta^2)*(1 - A/2/B);
    iV = eye(3) - 1/2*SO + Q*SO^2; % use directly skew of omega
end

omega = [SO(3,2) SO(1,3) SO(2,1)];
assert(all(isnan(omega) == 0),'not nan omega');                      
u = iV*t;

y = [omega,u(:)']';
