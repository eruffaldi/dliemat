function y = se3_fromRvecT(r,t)

y = se3_exprot(r);
y(1:3,4) = t;
