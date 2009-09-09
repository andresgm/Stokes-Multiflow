% Determina los nodos vecinos a cada nodo
function [n2nsparse,ia,ja] = node2nodesp(elements)

numnodes = max(elements(:));
numelements = size(elements,1);
elelist = repmat([1:1:numelements]',[1 3]);

ia = zeros(numelements*3,1);
ja = zeros(numelements*3,1);
contador = 0;
for i=1:numnodes
   nodesneigh = elements(elelist(elements(:) == i),:);
   nodesneigh(nodesneigh == i) = [];
   temp1 = sort(unique(nodesneigh(:)));
   numvecinos = size(temp1,1);
   contadorant = contador + 1;
   contador = (contador + numvecinos);
   ia(contadorant:1:contador) = i; %ones(numvecinos,1).*i;
   ja(contadorant:1:contador) = temp1;
end

n2nsparse = sparse(ia,ja,1);
