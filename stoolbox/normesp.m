% Retorna la norma euclideana de las filas de una matriz
% NormEsp_i=(sum(j)Xij^2)^1/2
% formato de X: dim(Q,3) siendo Q el numero de vectores

function normesp = normesp(x)

normesp = (sum(x.^2,2)).^0.5;
