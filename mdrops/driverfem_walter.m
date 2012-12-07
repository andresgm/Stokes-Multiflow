clear; clc; close all;

%% opciones de carga de archivos
% nombre de archivo a cargar y carpeta
nombreorigen = 'sph ref 3.mat';

geom.gssk = 1;
geom.csk = 1;

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

%% Geometria de referencia
geom.ref = geom.nodes;
geom.dsref = geom.ds;

paropt.tipo = 'extended';
[geom.curv,geom.normal,geom.Kg] = curvparaboloid(geom,paropt);
[lapbelmat,geom.Kg] = discretelaplacebeltrami(geom);
[geom.curv] = curvlb(geom,lapbelmat);

% figure(1)
% trimesh(geom.elements,geom.nodes(:,1),geom.nodes(:,2),geom.nodes(:,3));
% hold on
% quiver3(geom.nodes(:,1),geom.nodes(:,2),geom.nodes(:,3),...
%     geom.normal(:,1),geom.normal(:,2),geom.normal(:,3));
% axis equal; view(90,0); xlabel('x1'); ylabel('x2'); zlabel('x3'); colorbar;
% getframe; title('area');

%% Calculo de la fuerza elastica mediante metodo de los elementos finitos

gssk = geom.gssk;
csk = geom.csk;

Mglobal = zeros(3*numnodes,3*numnodes);
rglobal = zeros(3*numnodes,1);

[xg,wg] = quadrature;
[Abase,Abasectr,Metricg,Metricgctr,Jacmet,Jacmetctr] = isobase(geom);

alphadef = 0.1;
geom.nodes = geom.nodes + geom.normal*alphadef;

[abase,abasectr,metricg,metricgctr,jacmet,jacmetctr] = isobase(geom);

tab = ...
    tabskalak(Metricgctr,Jacmetctr,metricg,metricgctr,jacmet,gssk,csk);
chi = zeros(3,3); % chi11, chi22 y chi12 en tres coordenadas cartesianas.

for i = 1:numelements
    for k = 1:7
        Na = p1basis(xg(k,1),xg(k,2));
        dNa = p1dbasis(xg(k,1),xg(k,2));
        wa = wg(k);
        wajac = wa*sqrt(jacmet(i));
        for p = 1:3
            pglob = 3*(geom.elements(i,p)-1)+1;
            chi(1,:) = dNa(p,1).*abase(:,1,i);
            chi(2,:) = dNa(p,2).*abase(:,2,i);
            chi(3,:) = 0.5*(dNa(p,1).*abase(:,2,i)+dNa(p,2).*abase(:,1,i));
            rglobal(pglob) = rglobal(pglob) ...
                + wajac*(chi(1,1)*tab(i,1)+chi(2,1)*tab(i,2)+2*chi(3,1)*tab(i,3));
            rglobal(pglob+1) = rglobal(pglob+1) ...
                + wajac*(chi(1,2)*tab(i,1)+chi(2,2)*tab(i,2)+2*chi(3,2)*tab(i,3));
            rglobal(pglob+2) = rglobal(pglob+2) ...
                + wajac*(chi(1,3)*tab(i,1)+chi(2,3)*tab(i,2)+2*chi(3,3)*tab(i,3));
            Nap = Na(p);
            for q = 1:3
                qglob = 3*(geom.elements(i,q)-1)+1;
                mele = wajac*Nap*Na(q);
                Mglobal(pglob,qglob) = Mglobal(pglob,qglob) + mele;
                Mglobal(pglob+1,qglob+1) = Mglobal(pglob+1,qglob+1) + mele;
                Mglobal(pglob+2,qglob+2) = Mglobal(pglob+2,qglob+2) + mele;
            end
        end
    end
end

qloadv = Mglobal\rglobal;
qload = reshape(qloadv,[3,numnodes])';

% % Esfuerzos normales
% isonorm = repmat(sum(isotens'.*geom.normal,2),[1 3]).*geom.normal;
% % Esfuerzos de marangoni
% rdeltafmaran = isotens'-isonorm;
% 
% 
% %% Resultados
figure(2);
grafscfld(geom,normesp(qload));
hold on
quiver3(geom.nodes(:,1),geom.nodes(:,2),geom.nodes(:,3),...
    qload(:,1),qload(:,2),qload(:,3));
axis equal; view(90,0); xlabel('x1'); ylabel('x2'); zlabel('x3'); colorbar;
getframe; title('q - skalak');
% 
% figure(3);
% % trimesh...
% %(geom.elements,geom.nodes(:,1),geom.nodes(:,2),geom.nodes(:,3),isotens);
% % hold on
% % quiver3(geom.nodes(:,1),geom.nodes(:,2),geom.nodes(:,3),...
% %     rdeltafmaran(:,1),rdeltafmaran(:,2),rdeltafmaran(:,2));
% grafscfld(geom,normesp(rdeltafmaran));
% hold on
% quiver3(geom.nodes(:,1),geom.nodes(:,2),geom.nodes(:,3),...
%     rdeltafmaran(:,1),rdeltafmaran(:,2),rdeltafmaran(:,3));
% axis equal; view(90,0); xlabel('x1'); ylabel('x2'); zlabel('x3'); colorbar;
% getframe; title('tension en el plano');
% hold off
