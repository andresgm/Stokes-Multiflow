% Determina los nodos vecinos a cada nodo
function node2nodecell = node2node(elements)

numnodes = max(elements(:));
numelements = size(elements,1);
elelist = repmat([1:1:numelements]',[1 3]);

node2nodecell = cell(numnodes,1);
for i=1:numnodes
   nodesneigh = elements(elelist(elements(:) == i),:);
   nodesneigh(nodesneigh == i) = [];
   temp1 = sort(unique(nodesneigh(:)));
   node2nodecell{i} = temp1;
end
end