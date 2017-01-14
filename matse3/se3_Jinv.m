function [iA,JA] = se3_Jinv(A)

iA = se3_inv(A);

JA = [ dtranspose(3,3) zeros(3,9); kron(eye(3),-A(1:3,4)')    -A(1:3,1:3)'];

function Tmn = dtranspose(n,m)
            d = m*n;
            Tmn = zeros(d,d);
            
            i = 1:d;
            rI = 1+m.*(i-1)-(m*n-1).*floor((i-1)./n);
            I1s = sub2ind([d d],rI,1:d);
            Tmn(I1s) = 1;
            Tmn = Tmn';
        end