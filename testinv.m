%% INVERSION simulation
% Alternatively make a function from: SE3 -> position using unscented
A = zeros(6,6);
q = se3d_cov(T2);
A(1:3,1:3) = q(1:3,1:3);
%A = q; % FULL
AA = se3d_set([1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1],A);
%AA = se3d_set(se3d_mean(T2),A);
[xp,xi,wei] = se3d_sigmas(AA);
yp = zeros(size(xp));
% inverse!
for I=1:size(xp,3)
    yp(:,:,I) = inv(squeeze(xp(:,:,I)));
end
[iT,C] = se3d_unsigmas(yp,xi,wei);


iAA = se3d_inv(AA);
se3d_mean(iT)
se3d_mean(iAA)

se3d_cov(iT)
se3d_cov(iAA)

