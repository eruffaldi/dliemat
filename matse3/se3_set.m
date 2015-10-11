function y = se3_set(R,t)

y = eye(4);
y(1:3,1:3) = R;
y(1:3,4) = t;
