addpath /Users/eruffaldi/Dropbox/PaperSRGVIZ/bib/barfoot_tro14/
addpath matse3

q = zeros(6,1);
q(6) = 0.1;
q(2) = 0.1;
cholSigma = diag(q);
Sigma = cholSigma*cholSigma';

w = cholSigma * randn(6,1);
a = vec2tranwrap(w);
w = cholSigma * randn(6,1);
b = vec2tranwrap(w);


order = 1;
[bm,bS] = compound(a,Sigma,b,Sigma,order);
[em,eS] = compoundwrap(a,Sigma,b,Sigma,order);
bm-em
bS-eS
bS
eS
