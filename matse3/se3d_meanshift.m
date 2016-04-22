% estimate statistics from data: X -> SE3 dist
%
% This is the meanshift, if gk is 1 we obtain se3d_est
%
% Emanuele Ruffaldi 2016
% See Also: http://coewww.rutgers.edu/riul/research/papers/pdf/mspose.pdf
function y = se3d_meanshift(x,steps,gk,fx)

N = size(x,3);

if nargin < 3 || isempty(gk)
gk = x(:,:,1); % starting => random selection
end

for k=1:steps
	igk = se3_inv(gk);
	ma = zeros(1,6);
	w = 0;
	for i=1:N
		v_ik = se3_log(x(:,:,i)*igk);
		w_ik = fx(v_ik);
		w = w  + w_ik;
		ma = ma + w_ik*v_ik(:)';
	end
	ma = ma / w;
	gk = se3_exp(ma)*gk;
end

Sk = zeros(6,6);
igk = se3_inv(gk);
for i=1:N
	v_ik = se3_log(x(:,:,i)*igk);
	Sk = Sk + v_ik'*v_ik;
end
if N > 2
	Sk = Sk / (N-1); % unbiased estimator
end
y = se3d_set(gk,Sk);