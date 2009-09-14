% calcula el producto de un arreglo de m matrices [3x3] por un arreglo de m
% vectores fila.
% Mij*Vj (opt 1)
% Mij*Vi (opt 2)
% salida: m vectores fila.
function prodmatvect = matvect(matriz,vector,opt)

if nargin < 3
    opt = 1 % multiplicacion matriz vector columna
end


if opt == 1
    % Mij*Vj
    vectorl = repmat(permute(vector,[3 2 1]),[3 1 1]);
    prodmatvect = permute(sum(matriz.*vectorl,2),[3 1 2]);
elseif opt == 2
    % Mij*Vi
    vectorl = repmat(permute(vector',[1 3 2]),[1 3 1]);
    prodmatvect = permute(sum(matriz.*vectorl,1),[3 2 1]);
end
    
    