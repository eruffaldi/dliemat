function Y = vec2tranwrap(u)

% reverse
% Y = se3_exp([u(4:6);u(1:3)]);
%
% SAME
Y = se3_exp(u);
