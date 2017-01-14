function [AB,JA,JB] = se3_Jmul(A,B)

AB = se3_mul(A,B);
JA = kron(B',eye(3));
JB = kron(eye(4),A(1:3,1:3));
