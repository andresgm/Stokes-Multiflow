clear;clc; %close all;

%% opciones de carga de archivos
% nombre de archivo a cargar y carpeta
nombreorigen = 'sph ref 3.mat';

noiseint = 0.05;
noiserep = 4;
def = .01;

geom.ks = 1e6;
geom.mu = 1;

%% Carga de geometria y parametros de simulacion
sbar = systembar();
load([cd sbar nombreorigen]);
% PROCESAMIENTO DE LA MALLA ORIGINAL
% Volumen original
Radius = max(normesp(Nodes));
% Numero de Elementos y numero de Nodos
geom.nodes = Nodes;
geom.elements = Elements;
geom.numnodes = size(geom.nodes,1);
geom.numelements = size(Elements,1);
numnodes = geom.numnodes;
numelements = geom.numelements;

% calcule el volumen inicial de la gota
normalandgeoopt.normal = 1;
normalandgeoopt.areas = 1;
normalandgeoopt.vol = 1;
geomprop = normalandgeo(geom,normalandgeoopt);
% geom.normalele = geomprop.normalele;
geom.normal = geomprop.normal;
geom.dsi = geomprop.dsi;
geom.ds = geomprop.ds;
geom.s = geomprop.s;
geom.vol = geomprop.vol;
geom.jacmat = geomprop.jacmat;
geom.volini = geom.vol;
geom.areaini = geom.s;

geom.nodecon2node = node2node(geom.elements);
lmin = zeros(numnodes,1);
for k = 1:numnodes
   % extraiga los nodos vecinos a un nodo en la misma gota 
   nodesadj = geom.nodecon2node{k};
   lmin(k) = min(normesp(repmat(geom.nodes(k,:),[size(nodesadj,1) 1])...
      - geom.nodes(nodesadj,:)));  
end

%% Ruido
% Longitud minima para verificacion de paso de tiempo.
lmint = min(lmin);

for k = 1:noiserep
    noisevel = ones(size(geom.nodes))...
            .*(rand(size(geom.nodes))-0.5)*lmint*noiseint;
    noisenormal = repmat(sum(noisevel.*geom.normal,2),[1 3]).*geom.normal;
    noisetan = noisevel - noisenormal;

    geom.nodes = geom.nodes + noisetan;
end

geomprop = normalandgeo(geom,normalandgeoopt);
% geom.normalele = geomprop.normalele;
geom.normal = geomprop.normal;
geom.dsi = geomprop.dsi;
geom.ds = geomprop.ds;
geom.s = geomprop.s;
geom.vol = geomprop.vol;
geom.jacmat = geomprop.jacmat;
geom.volini = geom.vol;
geom.areaini = geom.s;

%% Geometria de referencia
geom.ref = geom.nodes;

geom.nodes = geom.nodes + def*geom.normal;

[lapbelmat,geom.Kg] = discretelaplacebeltrami(geom);
[geom.curv] = curvlb(geom,lapbelmat);

%% Calculo de la fuerza elastica mediante metodo de los elementos finitos
[geom.shapeA, geom.shapeB, geom.refrot] = shapefun(geom);
[isotens] = isotension(geom,geom.ks);


%% Resultados
figure(1);
grafscfld(geom,isotens);
axis equal; view(90,0); xlabel('x1'); ylabel('x2'); zlabel('x3'); colorbar;
getframe; title('area');
