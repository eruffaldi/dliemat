%   n     - Dimensionality of random variable
%   alpha - Transformation parameter  (optional, default 0.5)
%   beta  - Transformation parameter  (optional, default 2)
%   kappa - Transformation parameter  (optional, default 3-n)
% Copyright (C) 2006 Simo Sarkka
%
% Modified Emanele Ruffaldi 2014
function [wei] = ut_weights(n,alpha,beta,kappa)

lambda = alpha^2*(kappa+n)-n;

WM = repmat(1/(2*(n+lambda)), 2*n+1,1); % except first
WM(1) = lambda / (n+lambda); % first different

WC = WM;
WC(1) = lambda / (n+lambda) + (1-alpha^2+beta);

W = eye(length(WC)) - repmat(WM,1,length(WM));
W = W*diag(WC)*W';


wei = [];
wei.WC = WC;
wei.WM = WM;
wei.W = W;

