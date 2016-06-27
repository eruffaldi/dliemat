%
% Product/Composition
%
% Emanuele Ruffaldi SSSA 2015
%
% Fourth order Method from 
% From: Timothy D Barfoot and Paul T Furgale, 
%       Associating Uncertainty with Three-Dimensional Poses for use in Estimation Problems
%		DOI: 10.1109/TRO.2014.2298059
%       tim.barfoot@utoronto.ca, paul.furgale@mavt.ethz.ch

function y=se3d_mul(a,b,order)

if nargin == 2
    order = 2;
end

[ga,ca] = se3d_get(a);
[gb,cb] = se3d_get(b);
A = se3_adj(ga);
S2p = A*cb*A';

if order == 4
    
   % Fourth-order method from barfoot
   Sigma1rr = ca(1:3,1:3);
   Sigma1rp = ca(1:3,4:6);
   Sigma1pp = ca(4:6,4:6);
   
   Sigma2rr = S2p(1:3,1:3);
   Sigma2rp = S2p(1:3,4:6);
   Sigma2pp = S2p(4:6,4:6);
   
   A1 = [covop1(Sigma1pp) covop1(Sigma1rp+Sigma1rp'); zeros(3) covop1(Sigma1pp) ];
   A2 = [covop1(Sigma2pp) covop1(Sigma2rp+Sigma2rp'); zeros(3) covop1(Sigma2pp) ];
   
   Brr = covop2(Sigma1pp,Sigma2rr)+covop2(Sigma1rp',Sigma2rp)+covop2(Sigma1rp,Sigma2rp')+covop2(Sigma1rr,Sigma2pp);
   Brp = covop2(Sigma1pp,Sigma2rp')+covop2(Sigma1rp',Sigma2pp);
   Bpp = covop2(Sigma1pp,Sigma2pp);
   B = [Brr Brp; Brp' Bpp];
   
   S = ca + S2p + (A1*S2p + S2p*A1' + A2*ca + ca*A2')/12 + B/4;

else
   S = ca + S2p;

end

y = se3d_set(ga*gb,S); % instead of: se3_mul(ga,gb)


% covariance operation 1
function A = covop1( B )


A = -trace(B)*eye(3) + B;


% covariance operation 2
function A = covop2( B, C )


A = covop1(B)*covop1(C) + covop1(C*B);




