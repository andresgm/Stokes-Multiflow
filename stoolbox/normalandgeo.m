function geomprop = normalandgeo(geom,normalandgeoopt,opcion)
% calcula la normal promedio en los nodos de una malla triangular
% calcula el area superficial de cada elemento y total
% calcula el area baricentrica alrededor de cada nodo
% calcula el volumen total del cuerpo
%
% InPUts: (geom,normalandgeomopt)
% - geom: structura compuesta mInImO por
% geom.nodes: array de coordenadas dim(numnodes,3)
% geom.elements: array de elementos dim(numelement,3)
% geom.element2node: cellarray de conectividad de nodos con elementos
% ver funcion preprogeom()
% 
% - normalandgeoopt: structura de opt calculo 1: caLcULe 0: nO caLcULe
% normalandgeoopt.normal: calcule normal promedio en cada nodo de
% geom.nodes
% normalandgeoopt.areas: calcule el area superficial de cada elemento de
% geom.elements y el area baricentrica alrededor de cada nodo de geom.nodes
% ademas del areas superficial total del cuerpo.
% normalandgeoopt.vol: calcule el volumen del cuerpo (debe activar
% normalandgeoopt.normal = 1 & normalandgeoopt.areas = 1
% 
% OUtPUts: (geom)
% en la misma estructura de entrada geom se adicionan las propiedades
% calculadas: 
% - geom
% geom.normal: normal en cada nodo dim(numnodes,3)
% geom.ds: area superficial de cada elemento dim(numelements,1)
% geom.dsi: area baricentrica alrededor de cada nodo dim(numnodes,1)
% geom.s: area superficial total del enmallado
% geom.vol: volumen total del enmallado

% calcule la metrica de transformacion a cada punto
if nargin < 3
    jacomp = metrictrans(geom,[1/3;1/3]);
else
    jacomp = metrictrans(geom,[1/3;1/3],opcion);
end

geomprop.jacmat = jacomp.jacinv;
geomprop.g = jacomp.g;

if isfield(normalandgeoopt,'normal') == 0
   normalandgeoopt.normal = 0;
end

if isfield(normalandgeoopt,'areas') == 0
   normalandgeoopt.areas = 0;
end

if isfield(normalandgeoopt,'vol') == 0
   normalandgeoopt.vol = 0;
end

if normalandgeoopt.normal == 1
   [geomprop.normal,geomprop.normalele] = normal(geom,jacomp);
elseif normalandgeoopt.normal == 0 && isfield(geom,'normal') == 0
   warning('No hay normal calculada y no se solicito el calculo')
else
   geomprop.normal = geom.normal;
end

if normalandgeoopt.areas == 1 
    [geomprop.s,geomprop.ds,geomprop.dsi] = areas(geom);
end
   
if normalandgeoopt.vol == 1 && normalandgeoopt.areas == 0
    error('To calculate volume surface areas must be computed use normalandgeoopt.areas = 1');
elseif normalandgeoopt.vol == 1 && normalandgeoopt.areas == 1
    parfield.jacomp = jacomp;
    parfield.normal = geomprop.normal;
    parfield.dsi = geomprop.dsi;
    volumen = volume(geom,parfield);
    geomprop.vol = volumen;
end

