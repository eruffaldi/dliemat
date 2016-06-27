function [T,S] = compoundwrap(t1,s1,t2,s2,xorder) 

if xorder == 1
    order = 2;
else
    order = 4;
end

% inconvert the incoming covariance to our format and then convertback
[T,S] = se3d_get(se3d_mul(se3d_set(t1,flipcov(s1)),se3d_set(t2,flipcov(s2)),order));
S = flipcov(S);

% barfoot14 covariance scheme is based on translation+rotation 
% our covariance scheme is rotation+
function S = flipcov(S)

P = [0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1; 1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0];
S = P*S*P;
