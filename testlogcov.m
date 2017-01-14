% T ~ N(mu,Sigma) ==== xi ~ N(log(mu),Xigma)
%
% Demonstrate that: Sigma = J(mu) Xigma J(mu)' using Sigma points
%
% 2016-07-01
addpath /Users/eruffaldi/Dropbox/repos/liemat/matse3
addpath /Users/eruffaldi/Dropbox/PaperSRGVIZ/bib/barfoot_tro14/
%T = se3d_set(se3_exp([0.2,0.3,0.4,0.5,0.1,0.0]),0.5*diag([0.1,0.7,0.2,0.4,0.5,0.6]));
M = vec2tran([0.2,0,0.2,0,0.1,0.0]');
T = se3d_set(M,0.01*diag([0.1,0.0,0.0,0.0,0.0,0.0]));
se3d_log = @tran2vec;

% Compute it using Sigma points
[xp,xi,wei] = se3d_sigmas(T);
yp = zeros(6,size(xp,3));
for I=1:size(xp,3)
    yp(:,I) = se3d_log(xp(:,:,I));
end

ym = yp*wei.WM;  
yC = yp*wei.W*xi;
yS = yp*wei.W*yp';

J = vec2jacInv(tran2vec(M));
ySS = J*se3d_cov(T)*J';
ym-tran2vec(se3d_mean(T))  % 1E-16, mean is good now
yS
ySS
diffCov = yS-ySS   % 1E-16