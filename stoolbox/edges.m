% Determina los nodos que define cada edge de la malla elements
function edgeindex = edges(elements,ele2node)

if nargin < 2
    % determine los elementos vecinos a cada nodo
    ele2node = element2node(elements);
end

numnodes = max(elements(:));
edgeindex = zeros(0,2);
for i = 1:numnodes
   eleindextemp = ele2node{i};
   nodecon2node = elements(eleindextemp,:);
   nodecon2node(nodecon2node == i) = [];
   temp1 = sort(unique(nodecon2node(:)));
   de = [repmat(i,[size(temp1,1) 1]) temp1];
   edgeindex = [edgeindex; de];
end
