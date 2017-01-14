function [r,JA,jp] = se3_Jinvtx(A,p)

%( I3 kron ((p-t_A)') | -R_A') 
r = se3_tx(se3_inv(A),p); % optimize
JA = kron(eye(3), [ [p(:)-A(1:3,4)' , - A(1:3,1:3)']]);
Jp = A(1:3,1:3)';