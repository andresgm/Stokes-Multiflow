% calcula el area superficial de un arreglo triangular
function [s,ds,dsi] = areas(geom)

if isfield(geom,'element2node') ~= 1
    % no hay tabla de elementos-nodos calculela
    ele2node = element2node(geom.elements);
else
    ele2node = geom.element2node;
end

numnodes = size(geom.nodes,1);

nodo1 = geom.nodes(geom.elements(:,1),:);
nodo2 = geom.nodes(geom.elements(:,2),:);
nodo3 = geom.nodes(geom.elements(:,3),:);

ds = trianglearea(nodo1,nodo2,nodo3);
if isfield(geom,'numdrops') == 1
   if geom.numdrops == 1 
        s = sum(ds,1);
   else
        s = zeros(geom.numdrops,1);
        for j=1:geom.numdrops
           % extraiga los elementos de la j malla
           s(j) = sum(ds(geom.neledrop(j,1):geom.neledrop(j,2)),1);
        end
   end
else
    s = sum(ds,1);
end
% asigne la 3ra parte a cada nodo y sume (area baricentrica)
dsi = zeros(numnodes,1);
parfor i = 1:numnodes
   % extraiga los elementos del i nodo
   iele = ele2node{i};
   dsi(i) = (1/3)*sum(ds(iele),1);
end



