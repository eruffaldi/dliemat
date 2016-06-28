function y = se3d_fromRvec(r,avar)
%
% Rotation Vector (axis by angle) and variance of the angle
T = se3_fromRvecT(r,[0,0,0]);
qn = norm(r);
S = zeros(6);
if qn > 0
    a = abs(r)./norm(r); % FIX ME
    S(4:6,4:6) = avar*diag(a);
end
y = se3d_set(T,S);

