% proyecta de manera iterativa un elipsoide con volumen igual al de una
% esfera unitaria (4/3)pi y exceso de area dado por excesoareaobj y lo
% guarda en un archivo sph num2save.mat. el valor de num2save debe ser
% mayor que 5. La toleracia del error de exceso de area es de 1e-5

% generacion de elipsoides - semilla
clc;
clear;

reflevel = 3;
excesoareaobj = 1;
num2save = 10;

a = 1;
b = 1;
c = .97;

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
geom.volini = geomprop.vol;
geom.areaini = geomprop.s;

errorvoltol = 1e-6;
optesc.maxit = 15000;
optesc.kp = 20;
optesc.deltate = 0.01;
optesc.tolerrorvol = errorvoltol;
geom = escaling(geom,optesc,1,.1);
geom.vertices = extractvertices(geom);

for i = 1:geom.numnodes
    if geom.nodes(i,3) >= 0
        geom.nodes(i,3) = c*(1-geom.nodes(i,1).^2-geom.nodes(i,2).^2);
    else
        geom.nodes(i,3) = -c*(1-geom.nodes(i,1).^2-geom.nodes(i,2).^2);
    end
end

Nodes = geom.nodes;
Elements = geom.elements;

grafscfld(geom,1); xlabel('x1'); ylabel('x2'); zlabel('x3'); view(90,0); axis equal
nombrearchivo = ['sph ref ' num2str(num2save) '.mat'];
save(nombrearchivo,'Nodes','Elements');
disp(['archivo ' nombrearchivo ' guardado'])

