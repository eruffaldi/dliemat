%
% Sets the distribution specifying mean and variance
%
% Internally the group is stored  as flattened
function y = se3d_set(g,S)

assert(length(g) == 4,"group has to be in 4x4 form")
assert(length(S) == 6,"covariance has to be in 6x6 form")

y = [g(:);S(:)];