%
% Logarithm of the group: the algebra
%
% Emanuele Ruffaldi 2016
%
function y = se3_log(x)

R = x(1:3,1:3);
t = x(1:3,4);

% this branch is from barfoot the others ethan eade report
if 1==0

    [v,d]=eig(R);
    for i=1:3 
       if abs(d(i,i)-1) < 1e-10 
          a = v(:,i);
          a = a/sqrt(a'*a);
          phim = acos((trace(R)-1)/2);
          omega = phim*a;
          z = se3_exprot(omega);
          if abs(trace(z(1:3,1:3)'*R)-3) > 1e-14
             omega = -omega;
          end
       end
    end
    iV = vec2jacInv(omega);

else
   
    % first log of R
    %simplify(taylor(sin(x)/x,x,0,'Order',6))
    %       x^4/120 - x^2/6 + 1
    %simplify(taylor((1-cos(x))/x/x,x,0,'Order',6))
    %       x^4/720 - x^2/24 + 1/2
    %simplify(taylor((1-sin(x)/x)/x/x,x,0,'Order',6))
    %       x^4/5040 - x^2/120 + 1/6
    theta = acos((trace(R)-1)/2); %acos(max(-1,min((trace(R)-1)/2,1)));
    if isreal(theta) == 0
        R = R/abs(det(R));
        theta =  acos((trace(R)-1)/2);
    end

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


end

u = iV*t;
y = [omega(:);u(:)];


function invJ =vec2jacInv(phi)


    tolerance = 1e-12;
    ph = norm(phi);
    if ph < tolerance
        % If the angle is small, fall back on the series representation
        invJ = vec2jacInvSeries(phi,10);
    else
        axis = phi/norm(phi);
        ph_2 = 0.5*ph;

        invJ =   ph_2 * cot(ph_2)* eye(3)...
               + (1 - ph_2 * cot(ph_2))* axis * axis'...
               - ph_2 * hat(axis);
    end   


function S = hat(v)
S = [  0   -v(3)  v(2)
    v(3)  0    -v(1)
    -v(2) v(1)   0];


        
function [ invJ ] = vec2jacInvSeries( vec, N )
    
    invJ = eye(3);
    pxn = eye(3);
    px = hat(vec);
    for n = 1:N
        pxn = pxn * px/n;
        invJ = invJ + bernoullinumber(n) * pxn;
    end
    

