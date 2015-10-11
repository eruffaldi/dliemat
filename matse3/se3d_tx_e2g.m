% f is a function from euclidean to group
%
% input is covariance m,S in K dimensional space
% output is SE3 Gaussian
function [yd,xyS] = se3d_tx_e2g(m,S,f,sparams)

% TODO compute sigma points for m,S

% evaluate
yp = zeros(size(xp,1),4,4);
for I=1:size(xp,1)
	yp(I,:,:) = f(xp(I,:)));
end

[yd,xyS] = se3d_unsigmas(yp,X,wei);

