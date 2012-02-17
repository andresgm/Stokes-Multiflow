clear; clc; close all;

%% opciones de carga de archivos
% nombre de archivo a cargar y carpeta
nombreorigen = 'sph ref 3.mat';

noiseint = 0.025;
noiserep = 20;
def = -1e-2;

geom.ks = 1e6;
geom.mu = 0;

rkmaran = 1;

%% Carga de geometria y parametros de simulacion
sbar = systembar();
load([cd sbar nombreorigen]);
% PROCESAMIENTO DE LA MALLA ORIGINAL
% Numero de Elementos y numero de Nodos
geom.nodes = Nodes;%[0 0 0; 1 0 0; .5 .5 0 ; 0 1 0; 1 1 0];
geom.elements = Elements;%[1 2 3; 2 5 3; 5 4 3;4 1 3];
geom.numnodes = size(geom.nodes,1);
geom.numelements = size(geom.elements,1);
numnodes = geom.numnodes;
numelements = geom.numelements;

% calcule el volumen inicial de la gota
normalandgeoopt.normal = 1;
normalandgeoopt.areas = 1;
normalandgeoopt.vol = 1;
geomprop = normalandgeo(geom,normalandgeoopt,1);
geom.normalele = geomprop.normalele;%[0 0 1; 0 0 1;0 0 1;0 0 1]';%
geom.normal = geomprop.normal;%[0 0 1; 0 0 1;0 0 1; 0 0 1; 0 0 1];%
geom.dsi = geomprop.dsi;%[.2;.2;.2;.2;.2];%
geom.ds = geomprop.ds;%[.25;.25;.25;.25];%
geom.s = geomprop.s;%1;%
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

geomprop = normalandgeo(geom,normalandgeoopt,1);
geom.normalele = geomprop.normalele;
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
geom.dsref = geom.ds;
% 
[lapbelmat,geom.Kg] = discretelaplacebeltrami(geom);
[geom.curv] = curvlb(geom,lapbelmat);
paropt.tipo = 'extended';
[geom.curv,geom.normal,geom.Kg] = curvparaboloid(geom,paropt);

geom.nodes = geom.nodes + def*geom.normal;
% nombreorigen = 'ellipsoid.mat';
% sbar = systembar();
% load([cd sbar nombreorigen]);
% geom.nodes = Nodes;
% geom.nodes = [0 0 0; 1 0 0; .5 .5 .1 ; 0 1 0; 1 1 0];

[lapbelmat,geom.Kg] = discretelaplacebeltrami(geom);
[geom.curv] = curvlb(geom,lapbelmat);
paropt.tipo = 'extended';
[geom.curv,geom.normal,geom.Kg] = curvparaboloid(geom,paropt);

figure(1)
trimesh(geom.elements,geom.nodes(:,1),geom.nodes(:,2),geom.nodes(:,3));
hold on
quiver3(geom.nodes(:,1),geom.nodes(:,2),geom.nodes(:,3),...
    geom.normal(:,1),geom.normal(:,2),geom.normal(:,3));
axis equal; view(90,0); xlabel('x1'); ylabel('x2'); zlabel('x3'); colorbar;
getframe; title('area');

%% Calculo de la fuerza elastica mediante metodo de los elementos finitos
[geom.shapeA, geom.shapeB, geom.refrot] = shapefun(geom);
[isotens] = isotension(geom,geom.ks,geom.mu);
isonorm = repmat(sum(isotens'.*geom.normal,2),[1 3]).*geom.normal;

% Esfuerzos de marangoni
rdeltafmaran = isotens'-isonorm;


%% Resultados
figure(2);
grafscfld(geom,normesp(isonorm));
hold on
quiver3(geom.nodes(:,1),geom.nodes(:,2),geom.nodes(:,3),...
    isonorm(:,1),isonorm(:,2),isonorm(:,3));
axis equal; view(90,0); xlabel('x1'); ylabel('x2'); zlabel('x3'); colorbar;
getframe; title('tension normal');

figure(3);
% trimesh...
%(geom.elements,geom.nodes(:,1),geom.nodes(:,2),geom.nodes(:,3),isotens);
% hold on
% quiver3(geom.nodes(:,1),geom.nodes(:,2),geom.nodes(:,3),...
%     rdeltafmaran(:,1),rdeltafmaran(:,2),rdeltafmaran(:,2));
grafscfld(geom,normesp(rdeltafmaran));
hold on
quiver3(geom.nodes(:,1),geom.nodes(:,2),geom.nodes(:,3),...
    rdeltafmaran(:,1),rdeltafmaran(:,2),rdeltafmaran(:,3));
axis equal; view(90,0); xlabel('x1'); ylabel('x2'); zlabel('x3'); colorbar;
getframe; title('tension en el plano');
hold off
