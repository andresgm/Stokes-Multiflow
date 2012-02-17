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
s = sum(ds,1);
% asigne la 3ra parte a cada nodo y sume (area baricentrica)
% TODO: Revisar, esta si es el ?rea m?s conveniente? Voronoi?
dsi = zeros(numnodes,1);
for i = 1:numnodes
   % extraiga los elementos del i nodo
   iele = ele2node{i};
   dsi(i) = (1/3)*sum(ds(iele),1);
end



