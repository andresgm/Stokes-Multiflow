function dN = p1dbasis(eta1, eta2)
% Computes the value of the derivative of the basis for a certain value of
% intrinsic coordinates in an isoparametric triangular element

dN = zeros(3,2);

dN(1,:) = [-1,-1];
dN(2,:) = [1,0];
dN(3,:) = [0,1];

return