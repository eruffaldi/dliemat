function [r,JA,jp] = se3_Jtx(A,p)

r = se3_tx(A,p);
JA = kron([p(:)' 1],eye(3))
Jp = A(1:3,1:3);