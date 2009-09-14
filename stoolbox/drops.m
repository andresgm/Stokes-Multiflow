% Crea los nodos de las gotas adicionales a la gota establecida por defecto
% Modifica los nodos de tal manera que asigna los radios xr y desplaza la gota
% a las coordenadas xc
% establece un indice de los nodos y elementos que pertenecen a la j gota
function geom = drops(geom,xc,xr)
 
nodes = geom.nodes;
elements = geom.elements;
numelements = size(geom.elements,1);
numnodest = size(geom.nodes,1);
elemin = 1;
nodemin = 1;
geom.nodes= [];
geom.elements  = [];

geom.nnodesdrop = zeros(geom.numdrops,2);
geom.neledrop = zeros(geom.numdrops,2);
for j=1:geom.numdrops
   % extriga las coordendas del centroide de la j esima gota
   xcdrop = xc(j,:);
   %extraiga el radio de la j gota
   xrdrop = xr(j);
   % Mueva los nodos de la j gota a la posicion deseada (coordenadas)
   tempnodes = nodes*xrdrop + repmat(xcdrop,[size(nodes,1),1]);
   tempelements = elements + numnodest*(j-1);
   % Generar los nodos y los elementos
   geom.nodes = [geom.nodes;tempnodes];
   geom.elements = [geom.elements;tempelements];
   % Indice de nodos de la gota j 
   numnodes = size(geom.nodes,1);
   geom.nnodesdrop(j,:) = [nodemin numnodes];
   % Indice de elementos de la gota j
   numelements = size(geom.elements,1);
   geom.neledrop(j,:) = [elemin numelements];
   nodemin = numnodes + 1;
   elemin = numelements + 1;
end

geom.numnodes = numnodes;
geom.numelements = numelements;

