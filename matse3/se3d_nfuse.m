%
% Fuses N Gaussians using BarFoot14
%
% Parameters:
% x is the array of SE(3) distributions stacked by row
% steps is the number of iterations in the Gauss-Newton
% terms is the number of terms for the resolution of the inverse
% gk is the starting point
%
% Returns:
% 
%
function [y,V] = se3d_nfuse(x,steps,terms,gk)

if nargin < 4
    gk = T{1};
end

N = size(x,2);
invSigma = cell(N,1);
Sigma = cell(N,1);
cholSigma = cell(N,1);
T=cell(N,1);

% expand each in cella rrays
for k=1:N
    [R0,S0] = se3d_get(x(:,k));
    T{k} = R0;
    Sigma{k} = S0;
    invSigma{k} = inv( S0 );
    cholSigma{k} = chol( S0, 'lower' );
end



% for i=1:N
%     [Ri,Si] = se3d_get(x(:,i));
%     li = se3_log(Ri);
%     iSi = inv(Si);
%     iS = iS + iSi;
%     Q = Q + Si\li;
% end
% m0 = se3_exp(iS*Q); % eq 16
% S0 = inv(iS); % eq 28
%In these examples, the first-order constraint equations (21)
%and (22) were simultaneously solved by minimizing the sum of the squares of all of the
%elements on the left-hand side of these equations. If this minimization reaches a value of
%zero, we consider the constraints to have been solved.

   Test = gk;
   for i=1:steps      % Gauss-Newton iterations
      LHS = zeros(6);
      RHS = zeros(6,1);
      for k=1:N
         xik = se3_log( Test/T{k} );
         if terms <= 6
            invJ = vec2jacInvSeries( xik, terms );
         else
            invJ = vec2jacInv( xik );
         end
         invJtS = invJ'*invSigma{k};
         LHS = LHS + invJtS*invJ;
         RHS = RHS + invJtS*xik;
      end
      xi = -LHS \ RHS;
      assert(abs(det(Test(1:3,1:3))-1.0) < 1e-5,'a');
      Test = se3_exp( xi )*Test;
      assert(abs(det(Test(1:3,1:3))-1.0) < 1e-5,'b');
   end

   % How low did the objective function get?
   V = 0;
   for k=1:N
      xik = se3_log( Test*inv(T{k}) );
      V = V + xik'*invSigma{k}*xik/2;
   end

   % How close is the estimate to the true pose? NOT NEEDED
   %xi = se3_log( Ttrue*inv(Test) );
   %ep(n) = ep(n) + xi'*xi;

   % How big is the covariance? NOT
   Sigma_est = inv(LHS);
   %trcov(n) = trcov(n) + trace(Sigma_est'*Sigma_est);
   y = se3d_set(Test,Sigma_est);

end




function [ invJ ] = vec2jacInv( vec )
% VEC2JACINV Construction of the 3x3 J^-1 matrix or 6x6 J^-1 matrix in closed form
%   
% From: Timothy D Barfoot and Paul T Furgale, 
%       Associating Uncertainty with Three-Dimensional Poses for use in Estimation Problems
%		DOI: 10.1109/TRO.2014.2298059
%       tim.barfoot@utoronto.ca, paul.furgale@mavt.ethz.ch
%
%		Implementation of equation (99) or (103), depending on the input size, in the paper.
%
% input: 
%   vec: 3x1 vector or 6x1 vector
%
% output: 
%   invJ: 3x3 inv(J) matrix or 6x6 inv(J) matrix
%
   
%validateattributes(vec,{'double'},{'ncols',1})

tolerance = 1e-12;

if length(vec) == 3
    
    phi = vec;
    
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
    
elseif length(vec) == 6

    rho = vec(1:3);
    phi = vec(4:6);
    
    ph = norm(phi);
    if ph < tolerance;
        % If the angle is small, fall back on the series representation
        invJ = vec2jacInvSeries(vec,10);
    else
        invJsmall = vec2jacInv( phi );
        Q = vec2Q( vec );
        invJ = [ invJsmall -invJsmall*Q*invJsmall; zeros(3) invJsmall ];
    end
    
else   
    warning('vec2jacInv.m:  invalid input vector length\n');   
end   
   
end

function [ invJ ] = vec2jacInvSeries( vec, N )
% VEC2JACINVSERIES Construction of the 3x3 J^-1 matrix or 6x6 J^-1 matrix. Series representation
%   
% From: Timothy D Barfoot and Paul T Furgale, 
%       Associating Uncertainty with Three-Dimensional Poses for use in Estimation Problems
%		DOI: 10.1109/TRO.2014.2298059
%       tim.barfoot@utoronto.ca, paul.furgale@mavt.ethz.ch
%
%		Implementation of equation (99) or (103), depending on the input size, in the paper.
%
% input: 
%   vec: 3x1 vector or 6x1 vector
%   N:   number of terms to include in the series
%
% output: 
%   invJ: 3x3 inv(J) matrix or 6x6 inv(J) matrix
%

%validateattributes(vec,{'double'},{'ncols',1})
%validateattributes(N,{'numeric'},{'scalar'});

if length(vec)== 3 
    
    invJ = eye(3);
    pxn = eye(3);
    px = hat(vec);
    for n = 1:N
        pxn = pxn * px/n;
        invJ = invJ + bernoullinumber(n) * pxn;
    end
    
elseif length(vec) == 6    
    
    invJ = eye(6);
    pxn = eye(6);
    px = curlyhat(vec);
    for n = 1:N
        pxn = pxn * px/n;
        invJ = invJ + bernoullinumber(n) * pxn;
    end
    
else   
    warning('vec2jacInvSeries.m:  invalid input vector length\n');   
end

function [ b ] = bernoullinumber( k )
% BERNOULLINUMBER Generate the kth bernoulli number
% From: http://www.mathworks.com/matlabcentral/fileexchange/7257-bernoulli-numbers

if k == 0 
    b = 1;
elseif k == 1
    b = -1/2;
elseif mod(k,2)~= 0 
    b = 0;
else
    c = 1./factorial(2:k+1);
    b = (-1)^k .*factorial(k).* ...
    det(toeplitz(c, [c(1) 1 zeros(1, k-2)]));
end

end


end


function [ veccurlyhat ] = curlyhat( vec )
% CURLYHAT builds the 6x6 curly hat matrix from the 6x1 input
%
% From: Timothy D Barfoot and Paul T Furgale, 
%       Associating Uncertainty with Three-Dimensional Poses for use in Estimation Problems
%		DOI: 10.1109/TRO.2014.2298059
%       tim.barfoot@utoronto.ca, paul.furgale@mavt.ethz.ch
%
%		Implements equation (12) in the paper.  
%
% input: 
%   vec: 6x1 vector xi
%
% output: 
%   vechat: the 6x6 curly hat matrix 
%

%validateattributes(vec,{'double'},{'size',[6,1]});

phihat = hat( vec(4:6) );
veccurlyhat = [ phihat hat(vec(1:3)); zeros(3) phihat ];

end

function [ vechat ] = hat( vec )
% HAT builds the 3x3 skew symmetric matrix from the 3x1 input or 4x4 from 6x1 input
%
% From: Timothy D Barfoot and Paul T Furgale, 
%       Associating Uncertainty with Three-Dimensional Poses for use in Estimation Problems
%		DOI: 10.1109/TRO.2014.2298059
%       tim.barfoot@utoronto.ca, paul.furgale@mavt.ethz.ch
%
%		Implements equations (4) or (5), depending on input size, in the paper.
%
% input: 
%   vec: 3x1 vector phi or 6x1 vector xi
%
% output: 
%   vechat: the 3x3 skew symmetric matrix that can be used
%             to implement the cross product, or 4x4 for transformation
%             matrices
%

%validateattributes(vec,{'double'},{'ncols',1})

if size(vec,1) == 3
    
    vechat = [  0,     -vec(3),  vec(2);
            vec(3),   0    , -vec(1);
           -vec(2),  vec(1),   0    ];  
       
elseif size(vec,1) == 6
    vechat = [ hat( vec(4:6,1) ) vec(1:3,1); zeros(1,4) ];    
else   
    warning('hat.m:  invalid vector length for hat operator\n');   
end
end



function [ Q ] = vec2Q( vec )
% VEC2Q Construction of the 3x3 Q matrix
%   
% From: Timothy D Barfoot and Paul T Furgale, 
%       Associating Uncertainty with Three-Dimensional Poses for use in Estimation Problems
%		DOI: 10.1109/TRO.2014.2298059
%       tim.barfoot@utoronto.ca, paul.furgale@mavt.ethz.ch
%
%		Implementation of equation (102), in the paper.
%
% input:
%   vec: a 6x1 vector 
%
% output:
%   Q: the 3x3 Q matrix
%

%validateattributes(vec,{'double'},{'size',[6,1]});

rho = vec(1:3);
phi = vec(4:6);

ph = norm(phi);
ph2 = ph*ph;
ph3 = ph2*ph;
ph4 = ph3*ph;
ph5 = ph4*ph;

cph = cos(ph);
sph = sin(ph);

rx = hat(rho);
px = hat(phi);

t1 = 0.5 * rx;
t2 = ((ph - sph)/ph3) * (px*rx + rx*px + px*rx*px);
m3 = ( 1 - 0.5 * ph2 - cph ) / ph4;
t3 = -m3 * (px*px*rx + rx*px*px - 3*px*rx*px);
m4 = 0.5 * ( m3 - 3*(ph - sph - ph3/6)/ph5 );
t4 = -m4 * (px*rx*px*px + px*px*rx*px);

Q = t1 + t2 + t3 + t4;

end

