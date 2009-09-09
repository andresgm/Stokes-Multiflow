clear;clc
% system bar for path
sbar = systembar();

% driver de las rutinas
reflevel = 3;
filename = ['sph ref ' num2str(reflevel) '.mat'];
load([cd sbar filename]);
% Volumen original
Radius = max(normesp(Nodes));
% Numero de Elementos y numero de Nodos
geom.nodes = Nodes;
geom.elements = Elements;
numelements = size(geom.elements,1);
numnodes = size(geom.nodes,1);
%% pruebas individuales

ele2node = element2node(geom.elements);
[s,ds,dsi] = areas(geom);
% calcula los indices de los bordes
edgeindex = edges(geom.elements);
% nodos vecinos a cada nodo
nodecon2node = node2node(geom.elements);
% elementos vecinos a cada nodo
geom.element2node = element2node(geom.elements);
% elementos vecinos a cada elemento
ele2elecell = ele2ele(geom.elements);
% elementos vecinos a cada nodo -formato SPARSE-
geom.ele2nodessparse = ele2nodesp(geom.elements);

% metrica de transformacion 
jacomp = metrictrans(geom,[1/3;1/3]);
% calculo de la normal
[normal1,normalele1] = normal(geom);
% calculo de las areas
[s1,ds1,dsi1] = areas(geom);
% calculo del volumen
volumen1 = volume(geom);
% calculo de la curvatura con bestparaboloid
[curvatnodes,normalatnodes] = curvbestparaboloid(geom);
% calculo de la curvatura con laplace beltrami
tic
curvatnodes = curvlb(geom);
toc

[greeninffcn,closenode] = greeninf([0,0,0],geom.nodes,0);

% integrando del particle stress tensor
integrando1 = partstrtensor(randn(numnodes,3),geom.nodes,randn(numnodes,3),normal1,[1 0.1]);
% integral del trapecio para arreglos matriciales
intval = inttrapeciomata(dsi1,integrando1);

tic
[jacomp] = metrictrans(geom,[1/3;1/3]);
geom.element2node = element2node(geom.elements);
geom.nodecon2node = nodecon2node;
geom.edgeindex = edgeindex;
[normal2,normalele2] = normal(geom,jacomp);
[s2,ds2,dsi2] = areas(geom);


parfield.jacomp = jacomp;
parfield.normal = normal2;
parfield.dsi = dsi2;
volumen2 = volume(geom,parfield);
geom.normal = normal2;
tic
[curvatnodes,normalatnodes] = curvbestparaboloid(geom);
toc

e1 = s1 - s2;
e2 = ds1 - ds2;
e3 = dsi1 - dsi2;
e4 = volumen1 - volumen2;
e5 = normal1 - normal2;
e6 = normalele1 - normalele2;

% termino 1 de sulfactantes
velt = rand(size(geom.nodes,1),3);
tic
matterm1 = c_term1mat(geom,velt);
toc
% 

