%%
addpath purefx

%% Testing the group

g1 = se3_set(angle2dcm(pi/2,0,0),[1,0,0]);
g2 = se3_set(angle2dcm(0,pi/3,0),[0,1,0]);

se3_mul(g1,g2)-g1*g2 % zero
g1*se3_inv(g1) % identity 
se3_iadj(se3_adj(g1))-g1
se3_log(g1)


se3_log(se3_exp([1 0 1 0 0.5 0])) % wrong
se3_exp(se3_log(g1))-g1 % wrong

