function y = se3_adj(x)

R = x(1:3,1:3);
t = x(1:3,4);
y = zeros(6);
y(1:3,1:3) = R;
y(1:3,4:6) = skew(t)*R; 
y(4:6,4:6) = R;

function S = skew(v)
S = [  0   -v(3)  v(2)
      v(3)  0    -v(1)
     -v(2) v(1)   0];
