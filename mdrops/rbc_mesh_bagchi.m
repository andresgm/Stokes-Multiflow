% proyecta de manera iterativa un elipsoide con volumen igual al de una
% esfera unitaria (4/3)pi y exceso de area dado por excesoareaobj y lo
% guarda en un archivo sph num2save.mat. el valor de num2save debe ser
% mayor que 5. La toleracia del error de exceso de area es de 1e-5

% generacion de elipsoides - semilla
clc;
clear;

reflevel = 4;
name2save = 'rbc_ref4_ea4.1466';

c0 = 0.200;
c1 = 2.003;
c2 = -1.123;
R = 1.3841;

% cargue la esfera
fileload = ['sph ref ' num2str(reflevel) '.mat'];
% Cargue el archivo de solucion *.mat
load([cd '/' fileload]);

% numero de gotas
geom.numdrops = 1;
% Coordenadas de los centroides de las gotas
xc =[0 0 0];
% Introduzca el/los radios de la/s gotas
xr=1;
% numero de elementos y numero de nodos

geom.nodes = Nodes;
geom.elements = Elements;
geom.numnodes = size(geom.nodes,1);
geom.numelements = size(geom.elements,1);
geom.numnodes = size(geom.nodes,1);
numnodes = geom.numnodes;

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
geom.volini = geom.vol;
geom.areaini = geom.s;

% Calcular la coordenada z 

for i = 1:geom.numnodes
    eta = geom.nodes(i,1);
    zita = geom.nodes(i,2);
    r2 = (eta^2+zita^2);
    geom.nodes(i,1)=R*eta;
    geom.nodes(i,2)=R*zita;
    b = (c0 + c1*r2 + c2*(r2^2));
    Z = 0.5*R*(sqrt(abs(1.0-r2)))*b;
    
% Calcular la magnitud del vector

    x2 = geom.nodes(i,1)^2;
    y2 = geom.nodes(i,2)^2;
    z2 = geom.nodes(i,3)^2;

    mag = sqrt(x2+y2+z2);
    
  if ((geom.nodes(i,3)/mag)<0);
	geom.nodes(i,3)=-Z;
  else 
	geom.nodes(i,3)=Z;
  end

end

[lapbelmat,geom.Kg] = discretelaplacebeltrami(geom);
[geom.curv] = curvlb(geom,lapbelmat);
paropt.tipo = 'extended';
[geom.curv,geom.normal,geom.Kg] = curvparaboloid(geom,paropt);

normalandgeoopt.normal = 0;
normalandgeoopt.areas = 1;
normalandgeoopt.vol = 1;
geomprop = normalandgeo(geom,normalandgeoopt);
geom.dsi = geomprop.dsi;
geom.ds = geomprop.ds;
geom.s = geomprop.s;
geom.vol = geomprop.vol;


disp(['Initial Volumen: ',num2str(geom.volini)]);
disp(['Volumen: ',num2str(geom.vol)]);
errorvol = (geom.vol - geom.volini)./geom.volini;
disp(['Error in Volume: ',num2str(errorvol)]);
volred = 6*sqrt(pi)*geom.vol/geom.s^(3/2);
disp(['Volumen reducido: ',num2str(volred)]);
disp(['Area: ',num2str(geom.s)]);
R0 = (3*geom.vol/(4*pi))^(1/3);
exarea = geom.s/(R0^2) - 4*pi;
disp(['Excess Area: ',num2str(exarea)]);

Nodes = geom.nodes;
Elements = geom.elements;

grafscfld(geom,1); xlabel('x1'); ylabel('x2'); zlabel('x3'); view(90,0); 
axis equal;
nombrearchivo = [name2save '.mat'];
save(nombrearchivo,'Nodes','Elements');
disp(['archivo ' nombrearchivo ' guardado'])

