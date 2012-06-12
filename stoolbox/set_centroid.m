% Crea los nodos de las gotas adicionales a la gota establecida por defecto
% Modifica los nodos de tal manera que asigna los radios xr y desplaza la gota
% a las coordenadas xc
% establece un indice de los nodos y elementos que pertenecen a la j gota
function geom = set_centroid(geom,xc)
 
nodes = geom.nodes;

% Mueva los nodos a la posicion deseada (coordenadas)
geom.nodes = nodes + repmat(xc,[size(nodes,1),1]);


