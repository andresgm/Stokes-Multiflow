% calcula el vector normal exterior a los elementos y a los nodos
function [normnode,normele] = normal(geom,jacomp)

if nargin < 2
    % calcule la metrica de transformacion a cada punto
    jacomp = metrictrans(geom,[1/3;1/3]);
end

if isfield(geom,'element2node') ~= 1
    % no hay tabla de elementos-nodos calculela
    ele2node = element2node(geom.elements);
else
    ele2node = geom.element2node;
end

numnodes = size(geom.nodes,1);

% normal en cada elemento
temp1 = [jacomp.g1 jacomp.g2 jacomp.g3]; 
normele = temp1 ./ repmat(normesp(temp1),1,3);
% normal Promedio en cada nodo
normnode = zeros(numnodes,3);
for i = 1:numnodes
   % extraiga los elementos del i nodo
   iele = ele2node{i};
   % extraiga y promedie las normales de los elementos
   temp1 = sum(normele(iele,:),1);
   normnode(i,:) = temp1./normesp(temp1);   
end    
