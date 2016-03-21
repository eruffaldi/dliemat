%
% Samples the gaussian in SE3 for n points
%
% Emanuele Ruffaldi 2015-2015
function y = se3d_sample(x,n)

[x,S] = se3d_get(x);
cS = chol(S);

y = zeros(4,4,n);

for i=1:n
	y(4,4,i) = se3_exp(randn(6,6)*cS)*x;
end
