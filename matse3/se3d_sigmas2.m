% computes sigma points for the lie group using params
%
% g is the distribution in compact form
% params is an optional parameter with:
% - alpha
% - beta
% - kappa
%
% Output:
% xp are the 2N+1 by sizeof[group] sigmapoints, with N=6 and sizeof[group]=4,4
% v  are the increments
% wei are the weights as of ut_mweights
%
% Reassemble instruction (see barefoot)
%   OUTPUT AS SE3: use se3d_est
%   OUTPUT AS linear: easy
% Emanuele Ruffaldi SSSA 2015-2016

function [xp,v,wei] = se3d_sigmas(d1,d2,params)


n = 12; % dimensionality of the problem

if nargin == 2 || isempty(params)
	params = struct('alpha',0.5,'beta',2,'kappa',3-n);
end

[g1,S1] = se3d_get(d1);
[g2,S2] = se3d_get(d2);

S = blkdiag(S1,S2);
C = cholcov2(S); % 12x12 -> 6x6
k = size(C,1);

wei = ut_mweights2(n,k,params.alpha,params.beta,params.kappa);
c = wei.WC;


% first compute the sigma points stored ad 4x4 in this implementation
xp = zeros(2*k+1,4,4,2);
xp(1,:,:,1) = g1; % not weighted
xp(1,:,:,2) = g2; % not weighted
% for covariance
v = zeros(2*k+1,12,1);
for I=1:k
    psi = c(I)*C(I,:); % dimension 
	psi1 = psi(1:6);
    psi2 = psi(7:12);
	
    v(I+1,:) = psi;
    v(I+1+k,:) = -psi;
    
    xp(I+1,:,:,1) = se3_mul(se3_exp(psi1),g1); % weighted local motion
	xp(I+1+k,:,:,1) = se3_mul(-se3_exp(psi1),g1); % weighred local motion
	xp(I+1,:,:,2) = se3_mul(se3_exp(psi2),g2); % weighted local motion
	xp(I+1+k,:,:,2) = se3_mul(-se3_exp(psi2),g2); % weighred local motion
end

% Linear form
%X(1,:) = mu(:)';
%for I=1:L
%    X(I+1,:) = X(1,:) + c(I)*S(I,:);
%    X(I+1+L,:) = X(1,:) - c(I)*S(I,:);
%end
