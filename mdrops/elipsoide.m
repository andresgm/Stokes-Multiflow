% proyecta de manera iterativa un elipsoide con volumen igual al de una
% esfera unitaria (4/3)pi y exceso de area dado por excesoareaobj y lo
% guarda en un archivo sph num2save.mat. el valor de num2save debe ser
% mayor que 5. La toleracia del error de exceso de area es de 1e-5

% generacion de elipsoides - semilla
clc;
clear;

reflevel = 3
excesoareaobj = 1
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
geom.volini = geomprop.vol;
geom.areaini = geomprop.s;

errorvoltol = 1e-6;
optesc.maxit = 15000;
optesc.kp = 20;
optesc.deltate = 0.01;
optesc.tolerrorvol = errorvoltol;
geom = escaling(geom,optesc,1,.1);
geom.vertices = extractvertices(geom);
veladapt = zeros(geom.numnodes,3);

errorexcesoarea = (0 - excesoareaobj)/excesoareaobj;
tolarea = 1e-6;
kprop1 = 1e-2;
kprop2 = 5e-3;
kprop3 = 1e-2;
k = 0;
maxit = 100;
nodosori = geom.nodes;
eje1 = [1 0 0];
eje2 = [0 1 0];
eje3 = [0 0 1];

while abs(errorexcesoarea) > tolarea;
    k = k + 1
    % Calcule velocidad de desplazamiento
    velnode1 = repmat((geom.normal*eje1')*kprop1,[1 3])...
        .*geom.normal.*repmat((geom.normal*eje1'),[1 3]);
    velnode2 = -repmat((geom.normal*eje2')*kprop2,[1 3])...
        .*geom.normal.*repmat((geom.normal*eje2'),[1 3]);
    velnode3 = -repmat((geom.normal*eje3')*kprop3,[1 3])...
        .*geom.normal.*repmat((geom.normal*eje3'),[1 3]);
    velnode = velnode1 + velnode2 + velnode3;
    velnormal = repmat(sum(velnode1.*geom.normal,2),[1 3]).*geom.normal;
    veladapt = meshadaptgrad(geom,velnormal,veladapt);
    f1 = (velnormal + veladapt);
    geom.nodes = geom.nodes + errorexcesoarea*f1;
    grafscfld(geom,1);
    xlabel('x1'); ylabel('x2'); zlabel('x3'); view(90,0); axis equal
    getframe;
    % invoque rutina de escalaje para escalar el elipsoide a igual volumen
    geomprop = normalandgeo(geom,normalandgeoopt);
    geom.normalele = geomprop.normalele;
    geom.normal = geomprop.normal;
    geom.dsi = geomprop.dsi;
    geom.ds = geomprop.ds;
    geom.s = geomprop.s;
    geom.vol = geomprop.vol;
    geom.jacmat = geomprop.jacmat;
    errorvol = abs((geom.vol - geom.volini)./geom.volini);
    geom = escaling(geom,optesc,1,errorvol);
    
    areafinal = geom.s
    volumen = geom.vol
    excesoarea = areafinal - geom.areaini;
    errorexcesoarea = (excesoarea - excesoareaobj)/excesoareaobj
    disp(['exceso area ' num2str(excesoarea)]);
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

