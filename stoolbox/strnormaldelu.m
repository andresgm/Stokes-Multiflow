function intdoublelayer = strnormaldelu(strfcn,normalv,velatnode)

numxpoles = size(strfcn,1);
intdoublelayer = cell(numxpoles,1);

m = max(size(normalv));
numnodes = m;

% entrada tipo lista
normalvvect = normalv;
% arreglo de valores de la normalv para ejecutar producto en paralelo
temp1 = permute(reshape(normalvvect',1,3*m),[1 3 2]);
normarray = reshape(repmat(temp1,[3 3 1]),[3 3 3 m]);
    
    
for i = 1:numxpoles
    % extraiga los tensores de esfuerzo para el ipole
    strfcnm = strfcn{i};
    % extraiga la velocidad del ipole
    velpole = velatnode(i,:);
    % Realice U(x) - U(x0)
    deltavel = velatnode - repmat(velpole,[numnodes 1]);

    % ejecute tijk nk para cada punto de integracion
    strnormalva = sum(strfcnm.*normarray,3);
    
   % permute para que quede array de 3 x 3 x numquadpoints*numelements
    strnormalva = permute(strnormalva,[1 2 4 3]);
    
   % invoque matvect para ejecutar dij*vi
    intdoublelayert = matvect(strnormalva,deltavel,2);
    
    intdoublelayer{i} = intdoublelayert;
end
    
