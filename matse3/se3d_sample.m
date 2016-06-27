%
% Samples the gaussian in SE3 for n points
%
% Emanuele Ruffaldi 2015-2015
function y = se3d_sample(x,n)

[x,S] = se3d_get(x);
ccS = cholcov(S);
cS = ccS'*ccS;

y = zeros(4,4,n);
Q = randn(6,n);
for i=1:n
    y(:,:,i) = se3_exp(cS*Q(:,i))*x;
end
    