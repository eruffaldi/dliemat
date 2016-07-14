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
function [y,Sxy] = se3d_unsigmas(yp,w,wei)

N = size(yp,3); % number of sigma points
v = zeros(6,N); % xi of output

steps = 20;

gk = squeeze(yp(:,:,1));  % average
% estimate mean but weighted of WM
for k=1:steps
	igk = se3_inv(gk);
	for i=1:N
		% same as: se3_logdelta but with igk once
        v(:,i) = se3_log(squeeze(yp(:,:,i))*igk); 
	end
	ma = v*wei.WM;	% 6 (2k+1)  (2k+1) 1
	gk = se3_exp(ma)*gk;
end

% last run for the Sk
igk = se3_inv(gk);
%imux = se3_inv(se3d_mean(gx));
for i=1:N
	v(:,i) = se3_log(squeeze(yp(:,:,i))*igk); 
	%w(i,:) = se3_log(xp(:,:,i)*inv(gx)));  % use from input
end
Sk = v*wei.W*v; % covariance
y = se3d_set(gk,Sk);

if nargout > 1
    Sxy = w*wei.W*v; % cross XY
end

% Non Matricial Form
% mu = mu + WM(i) * Y(:,i);
% S = S + WC(i) * (Y(:,i) - mu) * (Y(:,i) - mu)';
% C = C + WC(i) * (X(1:size(M,1),i) - M) * (Y(:,i) - mu)';
%
% Matricial Form
% S = Y*W*Y'
% C = X*W*Y'
% m = X*Wm;
