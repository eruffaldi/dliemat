%%
addpath matse3

%% Example 2
a = se3_fromRvecT([0.5,0,0],[0.2,0.3,0.0]);
b = se3_fromRvecT([0,0,0],[0.8,0.3,0.2]);
ea = se3_log(a);
aea = se3_exp(ea);
c = se3_mul(a,b);
c2 = a*b;
ec2 = se3_log(c2)
ia = se3_inv(a);
w = zeros(4,4,2);
w(:,:,1) = a;
w(:,:,2) = b;
e = se3d_est(w,5);
%% Testing the group

g1 = se3_set(angle2dcm(pi/2,0,0),[1,0,0]);
g2 = se3_set(angle2dcm(0,pi/3,0),[0,1,0]);

disp('Test manip of group')
se3_mul(g1,g2)-g1*g2 % zero
g1*se3_inv(g1) % identity 
se3_iadj(se3_adj(g1))-g1
se3_log(g1)
se3_exp(se3_log(g1))-g1 % wrong

disp('Test log/exp')
se3_log(se3_exp([0 0 0 0 0.5 0]))-[0 0 0 0 0.5 0]' % correct
se3_log(se3_exp([1 0 1 0 0 0]))-[1 0 1 0 0 0]' % correct
se3_log(se3_exp([1 0 1 0 0.5 0]))-[1 0 1 0 0.5 0]' % wrong

se3_log(se3_exp(se3_log(g1)))-se3_log(g1) % wrong


