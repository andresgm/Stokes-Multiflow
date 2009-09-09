% calcula el producto de un arraglo de matrices 3 x 3 x s por 1 vector 1 x 3
function producto = arraydot1vector(array,vector,opt)

if nargin < 3
    % por defecto aij*bj (suma sobre columnas)
    opt = 1;
end

[m,n,s] = size(array);

if opt == 1
    % por defecto aij*bj (suma sobre columnas)
    varray = repmat(vector,[3 1 s]);
    producto = permute(sum(array.*varray,2),[3 1 2]);    
elseif opt == 2
    % por defecto aij*bi (suma sobre filas)
    varray = repmat(vector',[1 3 s]);
    producto = permute(sum(array.*varray,1),[3 2 1]);    
end
