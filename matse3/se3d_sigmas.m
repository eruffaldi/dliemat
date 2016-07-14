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
% wi are the weights as of ut_mweights
%
% Note: for removing the mean term (zero mean assumption) clean up the
% weights
%
% Reassemble instruction (see barefoot)
%   OUTPUT AS SE3: use se3d_est
%   OUTPUT AS linear: easy
% Emanuele Ruffaldi SSSA 2015-2016

function [xp,v,wei] = se3d_sigmas(d,params)


n = 6; % dimensionality of the problem


[g,S] = se3d_get(d);

C = cholcov(S); % 6x6 -> 
k = size(C,1);

if nargin == 1 || isempty(params)
	params = struct('alpha',0.5,'beta',2,'kappa',3-k);
end

wei = ut_mweights2(n,k,params.alpha,params.beta,params.kappa);

% if k < n
%     m = 2*k+1;
%     wei.WC = wei.WC(1:m);
%     wei.WM = wei.WM(1:m);
%     wei.W = wei.W(1:m,1:m);
% end

c = wei.c;

% first compute the sigma points stored ad 4x4 in this implementation
xp = zeros(4,4,2*k+1);
xp(:,:,1) = g; % not weighted
v = zeros(6,2*k+1);
se3_exp = @vec2tran;
for I=1:k
    psi = c*C(I,:)';
    v(:,I+1) = psi;
    v(:,I+1+k) = -psi;
	xp(:,:,I+1) = se3_exp(psi)*g; % weighted local motion
	xp(:,:,I+1+k) = se3_exp(-psi)*g; % weighred local motion
end
v = v';

% Linear form
%X(1,:) = mu(:)';
%for I=1:L
%    X(I+1,:) = X(1,:) + c(I)*S(I,:);
%    X(I+1+L,:) = X(1,:) - c(I)*S(I,:);
%end
