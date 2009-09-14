% Determina los elementos alrededor de un nodo
function ele2node = element2node(elements)

nelements = size(elements,1);
nnodes = max(elements(:));

ele2node = cell(nnodes,1);

for i=1:nelements
    for j=1:3
        ele2node{elements(i,j)}(end + 1) = i;
    end
end
