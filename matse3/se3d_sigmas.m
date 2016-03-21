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

function [xp,v,wei] = se3d_sigmas(d,params)


n = 6; % dimensionality of the problem

if nargin == 1 || isempty(params)
	sparams = struct('alpha',0.5,'beta',2,'kappa',3-n);
end

[g,S] = se3d_get(d);

wei = ut_mweights(n,params.alpha,params.beta,params.kappa);
c = wei.WC;

C = chol(S); % TODO customize this function because is delicate

% first compute the sigma points stored ad 4x4 in this implementation
xp = zeros(2*n+1,4,4);
v = zeros(2*n+1,6);
xp(1,:,:) = g; % not weighted
v(1,:) = [0,0,0,0,0,0];
for I=1:n
	ee = se3_exp(c(I)*S(I,:));
	vI(I+1,:) = ee;
	vI(I+1,:) = -ee;
	xp(I+1,:,:) = se3_mul(g,ee)); % weighted local motion
	xp(I+1+n,:,:) = se3_mul(g,-ee); % weighred local motion
end

% Linear form
%X(1,:) = mu(:)';
%for I=1:L
%    X(I+1,:) = X(1,:) + c(I)*S(I,:);
%    X(I+1+L,:) = X(1,:) - c(I)*S(I,:);
%end
