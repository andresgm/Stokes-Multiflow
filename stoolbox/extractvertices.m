% Rutina que extrae todos los vertices de una gota
% Los vertices son necesarios en la rutina de adaptacion de malla.
function vertices = extractvertices(geom)

numnodes = geom.numnodes;
neighbors = geom.nodecon2node;
vertices = [];

disp('Encontrando vertices de la malla...')
tic
for i = 1:numnodes
   for j= 1:size(neighbors{i},1)
      if neighbors{i}(j) > i
         vertices(end+1,:) = [i,neighbors{i}(j)];
      end
   end
end
toc
disp(['Numero de vertices: ',num2str(size(vertices,1))])
disp(['Numero de elementos: ',num2str(geom.numelements)])
disp(['Numero de nodos: ',num2str(numnodes)])
end