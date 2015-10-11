% Reconstructs the SE3 Gaussian from the transformed sigma points yp 
% and computes the covariance using the original sigma points
%
% Assumes input and output SE3
%
% Inputs:
% yp transformed sigma points 2*N+1 by 6 with N=6
% w are the sigma point increments (originally xp is used as xp-mux) 
% gx is the SE3 Gaussian of the x (for computing covariance)
% wei is the weight structure used for computing the weights
function [y,Sxy] = se3d_unsigmas(yp,w,gx,wei)

N = size(x,3);
gk = (:,:,1);
Sk = zeros(6,6);
v = zeros(N,6);
w = zeros(N,6);

for k=1:steps
	igk = se3_inv(gk);
	for i=1:N
		% same as: se3_logdelta but with igk once
		v(i,:) = se3_log(se3_mul(yp(:,:,i)*igk)); 
	end
	ma = v*wei.Wm;	
	gk = se3_exp(ma)*gk;
end

% last run for the Sk
igk = se3_inv(gk);
%imux = se3_inv(se3d_mean(gx));
for i=1:N
	v(i,:) = se3_log(se3_mul(yp(:,:,i)*igk)); 
	%w(i,:) = se3_log(se3_mul(xp(:,:,i)*imux)); 
end
Sk = v*wei.W*v';
Sxy = w*wei.W*v';
y = se3d_set(gk,Sk);

% Non Matricial Form
% mu = mu + WM(i) * Y(:,i);
% S = S + WC(i) * (Y(:,i) - mu) * (Y(:,i) - mu)';
% C = C + WC(i) * (X(1:size(M,1),i) - M) * (Y(:,i) - mu)';
%
% Matricial Form
% S = Y*W*Y'
% C = X*W*Y'
% m = X*Wm;
