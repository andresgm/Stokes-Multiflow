% proyecta de manera iterativa un elipsoide con volumen igual al de una
% esfera unitaria (4/3)pi y exceso de area dado por excesoareaobj y lo
% guarda en un archivo sph num2save.mat. el valor de num2save debe ser
% mayor que 5. La toleracia del error de exceso de area es de 1e-5

% generacion de elipsoides - semilla
clc;
clear;

reflevel = 3
excesoareaobj = 0.3465
num2save = 10

if num2save < 6
    error('el numero del archivo a guardar debe ser mayor que 5')
end
% cargue la esfera
fileload = ['sph ref ' num2str(reflevel) '.mat'];
% Cargue el archivo de solucion *.mat
load([cd '/' fileload]);

% numero de gotas
geom.numdrops = 1;
% Coordenadas de los centroides de las gotas
xc =[0 0 0];
% Introduzca el/los radios de la/s gotas
xr=[1];
% numero de elementos y numero de nodos
geom.nodes = Nodes;
geom.elements = Elements;
geom.numnodes = size(geom.nodes,1);
geom.numelements = size(geom.elements,1);

geom = drops(geom,xc,xr);

% calcule las propiedades y tablas de la malla
% Index table
geom.indextable = 1:1:geom.numnodes;
% Tabla de elementos singulares a cada nodo
geom.element2node = element2node(geom.elements);
% Tabla de conectividad de nodos, bordes, e.t.c
geom.nodecon2node = node2node(geom.elements);
% Calcule las propiedades iniciales
normalandgeoopt.normal = 1;
normalandgeoopt.areas = 1;
normalandgeoopt.vol = 1;
geomprop = normalandgeo(geom,normalandgeoopt);
geom.normalele = geomprop.normalele;
geom.normal = geomprop.normal;
geom.dsi = geomprop.dsi;
geom.ds = geomprop.ds;
geom.s = geomprop.s;
geom.vol = geomprop.vol;
geom.jacmat = geomprop.jacmat;
geom.volini = 4.*pi/3.;
geom.areaini = 4.*pi;

errorvoltol = 1e-6;
optesc.maxit = 1000;
optesc.kp = 20;
optesc.deltate = 0.01;
optesc.tolerrorvol = errorvoltol;
geom = escaling(geom,optesc,1);

errorexcesoarea = 0 - excesoareaobj;
tolarea = 1e-6;
kprop = 0.8;
k = 0;
maxit = 1000;
nodosori = geom.nodes;
a = 1.0 - errorexcesoarea*kprop
b = 1.0 - errorexcesoarea*kprop
c = 1.0;

while abs(errorexcesoarea) > tolarea;
    k = k + 1;     
    % extrapole los nodos hacia elk elipsoide
    geom.nodes = [geom.nodes(:,1).*a geom.nodes(:,2).*b geom.nodes(:,3).*c];
    grafscfld(geom,1);
    xlabel('x1'); ylabel('x2'); zlabel('x3'); view(90,0); axis equal
    getframe;
    % invoque rutina de escalaje para escalar el elipsoide a igual volumen
    geom = escaling(geom,optesc,1);

    areafinal = geom.s
    volumen = geom.vol
    excesoarea = (areafinal - geom.areaini)/geom.areaini;
    errorexcesoarea = excesoarea - excesoareaobj
    disp(['exceso area ' num2str(excesoarea)]);
    a = 1.0 - errorexcesoarea*kprop
    b = 1.0 - errorexcesoarea*kprop
    if k == maxit
        error ('El escalaje no convergio')
    end
end

Nodes = geom.nodes;
Elements = geom.elements;

grafscfld(geom,1); xlabel('x1'); ylabel('x2'); zlabel('x3'); view(90,0); axis equal
nombrearchivo = ['sph ref ' num2str(num2save) '.mat'];
save(nombrearchivo,'Nodes','Elements');
disp(['archivo ' nombrearchivo ' guardado'])

