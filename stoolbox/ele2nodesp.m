% genere la matriz de conectividad de nodos vecinos en forma de matriz
% dispersa. 
% apropiada para hacer calculos donde intervienen nodos vecinos

function e2nsparse = ele2nodesp(elements)
% numero de nodos
numnodes = max(elements(:));

% genere el cell array de elementos vecinos
ele2nodescell = element2node(elements);
j = [];
i = [];
for k=1:numnodes
    % extraiga los elementos singulares al inode
    j = [j ele2nodescell{k}];
    i = [i ones(1,size(ele2nodescell{k},2)).*k];
end
e2nsparse = sparse(i,j,1);
