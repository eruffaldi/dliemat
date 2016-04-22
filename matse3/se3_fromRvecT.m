function y = se3_fromRvecT(r,t)
w = [r(:);0;0;0]';
y = se3_exp(w);
assert(abs(det(y(1:3,1:3))-1.0) <1e-5,'a');
y(1:3,4) = t;

