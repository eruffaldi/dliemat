% estimate statistics from data
function y = se3d_est(x,steps)

N = size(x,3);

gk = x(:,:,1);
Sk = zeros(6,6);

for k=1:steps
	igk = se3_inv(gk);
	ma = zeros(6,1);
	for i=1:N
		v_ik = se3_log(x(:,:,i)*igk);
		ma = ma + v_ik;
	end
	ma = ma / N;
	gk = se3_exp(ma)*gk;

end

Sk = zeros(6,6);
igk = se3_inv(gk);
for i=1:N
	v_ik = se3_log(x(:,:,i)*igk);
	Sk = Sk + v_ik*v_ik';
end
Sk = Sk / (N-1); % unbiased estimator
y = se3d_set(gk,Sk);