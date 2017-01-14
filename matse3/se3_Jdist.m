function [r,JA,JB] = se3_Jdist(A,B)

[iB,JiB] = se3_Jinv(B);
[AiB,JAiB_A,JAiB_iB] = se3_Jmul(A,iB);
[r,Jr] = se3_Jlog(AiB);

JA = Jr*JAiB_A;
JB = Jr*JAiB_iB*JiB;
