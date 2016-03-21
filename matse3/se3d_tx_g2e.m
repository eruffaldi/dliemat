% f is a function from the group to euclidean
%
% input is SE3 Gaussian, function and dimension of output
% output is euclidean in K dimensional space with covariances
%
% Emanuele Ruffaldi SSSA 2015-2016
function [ymu,yS,xyS] = se3d_tx_g2e(d,f,K,sparams)

[g,S] = se3d_get(d);

[xp,v,wei] = se3d_sigmas(g,S,sparams);

% then evaluate the function
yp = zeros(size(xp,1),K);
for I=1:size(xp,1)
	yp(I,:) = f(squeeze(xp(I,:,:)));
end

ymu = yp*wei.WM';
yS = yp*wei.W*yp.wei.W';
xyS = v*wei.W*Y'; 
