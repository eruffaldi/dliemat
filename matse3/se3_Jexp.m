function [A,J] = se3_Jexp(v)

% assume t omega
A = se3_exp(v);
J = [zeros(3) -skew([1 0 0]); zeros(3) -skew([0 1 0]); zeros(3) -skew([ 0 0 1]); eye(3) zeros(3)];

