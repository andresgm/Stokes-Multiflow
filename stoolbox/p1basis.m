function N = p1basis(eta1, eta2)
% Computes the value of the basis for a certain value of the intrinsic
% coordinates in an isoparametric triangular element

N = zeros(3,1);

N(1) = 1-eta1-eta2;
N(2) = eta1;
N(3) = eta2;

return