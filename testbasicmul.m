
%% SIMPLE compose simulation
% Alternatively make a function from: SE3 -> position using unscented
Q = se3_exp([0.2 0.3,0.4 0.1,0.2,0.3]);
Q = se3_exp([0 0 0 0 ,0.0,0.3]);
A = zeros(6,6);
q = se3d_cov(T2);
%A(1:3,1:3) = [1 0 0; 0 0 0; 0 0 0]; % q(1:3,1:3);
A(6,6) = 1.0;
A(5,5) = 1.0;
%A = q; % FULL
AA = se3d_set([1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1],A);
%AA = se3d_set(se3d_mean(T2),A);
[xp,xi,wei] = se3d_sigmas(AA);
yp = zeros(size(xp));
% inverse!
for I=1:size(xp,3)
    yp(:,:,I) = Q*squeeze(xp(:,:,I));
end
[iT,C] = se3d_unsigmas(yp,xi,wei);


iAA = se3d_mul(se3d_set(Q,zeros(6)),AA);
se3d_mean(iT)
se3d_mean(iAA)

se3d_cov(iT)
se3d_cov(iAA)
