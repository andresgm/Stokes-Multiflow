% Calcula un conjunto de propiedades geometricas de una malla superficial
% 3D
% 
% INPUTS:
% - geom: structura compuesta MINIMO por
% geom.nodes: array de coordenadas dim(numnodes,3)
% geom.elements: array de elementos dim(numelement,3)
% geom.element2node: cellarray de conectividad de nodos con elementos
% ver funcion preprogeom()
% geom.nodecon2node: cell array de conectividad de nodos con nodos
% ver funcion preprogeom()
% geom.edgeindex: array de conectividad de nodos con bordes
% ver funcion preprogeom()
%
% - propgeomopt: structura de opt calculo 1: CALCULE 0: NO CALCULE
% propgeomopt.normal: calcule la normal a cada nodo
%
% propgeomopt.curv: calcule la curvatura media en cada nodo
%   propgeomopt.curv = 'best parabolloid fitting' 
%       calcule la curvatura usando bestparabolloid fitting
%   propgeomopt.curv = 'parabolloid fitting' 
%       calcule la curvatura usando parabolloid fitting
%   propgeomopt.curv = 'laplace belt based' 
%       calcule la curvatura basada segun el op laplace beltrami
% propgeomopt.areas: calcule areas asociadas (Area de cada elemento, area
% baricentrica a cada nodo y area superficial)
%
% propgeomopt.vol: calcule volumen total del cuerpo
%
% propgeomopt.lapcurv.type: calcule el operador laplace beltrami sobre cada
% nodo en la superficie
% propgeomopt.lapcurv.type = 'voronoi' use area circuntentrica unicamente
% propgeomopt.lapcurv.type = 'bary' use area baricentrica unicamente
% propgeomopt.lapcurv.type = 'mixed' use los dos tipos de area
% (recomendado)
% propgeomopt.lapcurv.type = '' NO CALCULE LAPLACE BELTRAMI
%
% OUTPUTS:
% En la misma estructura de entrada geom se adicionan las propiedades
% calculadas: 
% - geom
% geom.normal: normal en cada nodo dim(numnodes,3)
% geom.ds: area superficial de cada elemento dim(numelements,1)
% geom.dsi: area baricentrica alrededor de cada nodo dim(numnodes,1)
% geom.s: area superficial total del enmallado
% geom.vol: volumen total del enmallado
% geom.curv: curvatura media en cada nodo
% geom.lapcurv: laplace beltrami de la curvatura media en cada nodo

function geom = propgeom(geom,propgeomopt)

% Recupere las opciones y variables

% Nodes = geom.nodes;
% Elements = geom.elements;

% Defina las opciones para calculo de la normal y areas
normalandgeoopt.normal = 0; 
normalandgeoopt.areas = 0; 
normalandgeoopt.vol = 0;
if propgeomopt.normal == 1
   normalandgeoopt.normal = 1; 
end
if propgeomopt.areas == 1 
   normalandgeoopt.areas = 1;
end
if propgeomopt.vol == 1 
   normalandgeoopt.vol = 1;
end

geom = normalandgeo(geom,normalandgeoopt);
    
if strcmp(propgeomopt.curv,'best parabolloid fitting') == 1
    % calcule la curvatura mediante best parabolloid fitting
geom = bestparaboloid(geom);   
% geom.curv = paraboloidcurv(geom);

% geom.curv = curvlb(geom);

elseif strcmp(propgeomopt.curv,'paraboloid fitting') == 1
    % calcule la curvatura mediante parabolloid fitting
    
elseif strcmp(propgeomopt.curv,'laplace belt based') == 1
    % calcule la laplace beltrami
    
end

if strcmp(propgeomopt.lapcurv.type,'voronoi') == 1
    lapbeltype = 'voronoi';
elseif strcmp(propgeomopt.lapcurv.type,'bary') == 1
    lapbeltype = 'bary';
elseif strcmp(propgeomopt.lapcurv.type,'mixed') == 1
    lapbeltype = 'mixed';    
else
    lapbeltype = '';
end

if strcmp(lapbeltype,'') == 0
    geom.lapcurv = laplacebeltrami(geom,geom.curv,lapbeltype);
else
    geom.lapcurv = 0;
end


