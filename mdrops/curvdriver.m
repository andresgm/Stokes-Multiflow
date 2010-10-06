clc; clear;

nombreorigen = 'sph ref 3';
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
geom.numdrops = 1;
xc =[0 0 0];
xr=[1];
geom = drops(geom,xc,xr);

geom.element2node = element2node(geom.elements);
geom.nodecon2node = node2node(geom.elements);
normalandgeoopt.normal = 1;
normalandgeoopt.areas = 1;
normalandgeoopt.vol = 1;
geomprop = normalandgeo(geom,normalandgeoopt);
geom.normalele = geomprop.normalele;
geom.normal = geomprop.normal;
clear Elements Nodes nombreorigen;

paropt.tipo = 'extended';
[geom.curv,geom.normal,geom.Kg] = curvparaboloid(geom,paropt);
geom.lapcurv = lapbel(geom,geom.curv);

lapbelmat = laplacebeltramimat(geom);
% geom.lapcurv2 = lapbelmat*geom.curv;
[geom.curv2] = curvlb(geom,lapbelmat);
geom.lapcurv3 = lapbelmat*geom.curv2;

figure(1);
% hold on
grafscfld(geom,geom.lapcurv);
axis equal; view(90,0); xlabel('x1'); ylabel('x2'); zlabel('x3'); colorbar;
getframe; title('laplace curv');
% 
% figure(2)
% grafscfld(geom,geom.lapcurv2);
% axis equal; view(90,0); xlabel('x1'); ylabel('x2'); zlabel('x3'); colorbar;
% getframe; title('laplace curv2');
% %quiver3(...
% %   geom.nodes(:,1),geom.nodes(:,2),geom.nodes(:,3),...
% %   geom.normal(:,1),geom.normal(:,2),geom.normal(:,3))
% %hold off
% 
figure(3);
grafscfld(geom,geom.curv);
axis equal; view(90,0); xlabel('x1'); ylabel('x2'); zlabel('x3'); colorbar;
getframe; title('Curv');
% 
% figure(4);
% grafscfld(geom,geom.Kg);
% axis equal; view(90,0); xlabel('x1'); ylabel('x2'); zlabel('x3'); colorbar;
% getframe; title('Gaussiana');

figure(1);
grafscfld(geom,geom.curv2);
axis equal; view(90,0); xlabel('x1'); ylabel('x2'); zlabel('x3'); colorbar;
getframe; title('Curv2');

figure(2)
grafscfld(geom,geom.lapcurv3);
axis equal; view(90,0); xlabel('x1'); ylabel('x2'); zlabel('x3'); colorbar;
getframe; title('laplace curv3');