function y = se3d_fuse(a,b)

[R0,S0] = se3d_get(a);
[R1,S1] = se3d_get(b);

Sy = S0 - S0*inv(S0 + S1)*S0;
v = se3_log(R1*se3_inv(R0));
A = Sy*inv(S1);
Ry = se3_exp(Sy*v)*R0; 

y = se3_set(Ry,Sy);

% OPTIMIZE ABOVE/VERIFY?
% v = R1*inv(R0)
% A = se3gadj(Sc*inv(S1))
% Rc = A*v*inv(A)*R0
