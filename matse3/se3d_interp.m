% interpolate 
% equivalent to slerp in pure quaternions
function y = se3d_interp(d1,d2,alpha)

[x1,S1] = se3d_get(d1);
[x2,S2] = se3d_get(d2);

ym = se3_interp(x1,x2,alpha);
yS = S1*alpha+(1-alpha)*S2;

y = se3d_set(ym,yS);