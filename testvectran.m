% Compare vectran implementation
addpath /Users/eruffaldi/Dropbox/PaperSRGVIZ/bib/barfoot_tro14/
addpath matse3

q = zeros(6,1);
q(6) = 0.1;
q(2) = 0.1;
cholSigma = diag(q);

w = cholSigma * randn(6,1);
a = vec2tran(w);
b = vec2tranwrap(w);
a-b
