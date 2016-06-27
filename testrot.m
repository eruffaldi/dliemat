% Build an SE3 distribution from a rotation around an axis angle with noise
% affected rotation

tn=2.0;
d = se3d_fromRvec(tn*[1 0.2 0],0.1); % build axis angle with variance along rotation
QQ = se3d_sample(d,1000); % sample the distribution

% extract
oo = zeros(size(QQ,3),3);
oon = oo;
oom = zeros(size(QQ,3),1);
for I=1:size(oo,1)
    l = se3_log(QQ(:,:,I)); % logarithm of transform
    oo(I,:) = l(1:3); % as rotation vector
    oom(I) = norm(oo(I,:)); % norm of the rotation
    oon(I,:) = oo(I,:)./oom(I); % normalized rotation vector
end

% back to distrbution via estimation
de = se3d_est(QQ,10);

% analyze averaging and norm
oon
mean(oom)
var(oom)

dem = se3d_mean(de)
ldem = se3_log(dem)
deS = se3d_cov(de)

%% Single Joint

d = se3d_fromRvec(2*[0 0 1],0.1); % Z rotation uncertain
d2 = se3d_set(se3_fromRvecT(eye(3),[1.0,0.0,0.0]),zeros(6)); % translate along x
de = se3d_mul(d,d2);
se3d_mean(de)
se3d_cov(de)

%%
d = se3d_fromRvec(tn*[1 0.2 0.3],0.1); % build axis angle with variance along rotation
Q = se3_fromRvecT([0.2,0.2,0.3],[1,2,3]);
dQ = se3d_set(Q,zeros(6));
dQ1 = se3d_mul(dQ,d);
dQ2 = se3d_mul(d,dQ);
disp('Resulting')
se3d_cov(dQ1)
se3d_cov(dQ2)
