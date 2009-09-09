% CALCULO DEL FLUJO DE STOKES PARA 2 GOTAS
% IMPLEMENTADO FLUJO INFINITO, SEMIINFINITO
% IMPLEMENTADO SINGLE LAYER

clear;clc
% Algoritmo de flujo de stokes con sulfactantes.
ca = 0.25;
lamda = 1;
% Adimensionalizacion
adim = 1;
% frecuencia de guardar resultados
outputfreq = 5;
% Banderas de fuerza dif 0: si. 1: no
ka = 1;
kb = 0;
kc = 0;
% malla a usar
reflevel = 2;
% numero de gotas
geom.numdrops = 2;
% Coordenadas de los centroides de las gotas
xc =[0 -1.2 0.25 ; 0 1.2 -0.25];
% xc =[0 0 0];
% Introduzca el/los radios de la/s gotas
xr=[1 1];
% xr=[1];

% pasos de tiempo de la simulacion
numtimesteps = 20000;
deltat = 0.001;
redfactor = 2;

% parametros de adaptacion
% velopt: 1 hidrodinamica velopt:2 normal
velopt = 1;
% meshadapt lowenberg
adaptparms.psi = 1;
adaptparms.lamda = lamda;
adaptparms.flow = flow;
% escalaje
errorvoltol = 1e-3;
optesc.maxit = 15000;
optesc.kp = 20;
optesc.deltate = 0.01;
optesc.tolerrorvol = errorvoltol;
% TODO: DEFINIR OPCIONES DE ESCALAJE
%% procesamiento de parametros
if adim == 1
    % adimensionalizacion de bazhlekov
    parms.rkextf = 2/(lamda + 1);
    parms.rksl = 2/(lamda + 1);
    parms.rkdl = 2*(lamda - 1)/(lamda + 1);

    % parametros de simulacion
    if ka == 1 
        % curvatura constante
        parms.rkcurv = 1/ca;
    else
        parms.rkcurv = 0;
    end

    if kb == 1
        % gravedad
        parms.rkgrav = 1/ca;
    else
        parms.rkgrav = 0;
    end    
elseif adim == 2
    
end

%% procesamiento de la malla
sbar = systembar();
% cargue el archivo base
filename = ['sph ref ' num2str(reflevel) '.mat'];
load([cd sbar filename]);
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
geom = drops(geom,xc,xr);

% Index table
% geom.indextable = [1:1:geom.numnodes];
    % Tabla de elementos singulares a cada nodo
geom.element2node = element2node(geom.elements);
%     % Elementos singulares a cada nodo tabla SPARSE
% geom.e2nsparse = ele2nodesp(geom.elements);
% TODO: Borrar de la rutina ele2nodesp
    % Tabla de conectividad de nodos, bordes, e.t.c
geom.nodecon2node = node2node(geom.elements);
% geom.edgeindex = edges(geom.elements);
% TODO: si no se necesita borrar edges
    
% calcule el volumen inicial de la gota
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
geom.volini = geom.vol;


%% evolucion mediante runge kutta de 2to orden
counter = 0;
geom.tiempo = 0;
itsaved = 0;
for p = 1:numtimesteps
% calcule la distancia minima de adaptacion y el paso de tiempo
    counter = counter + 1;    
   
    if geom.numdrops == 1
        % lmin entre nodos de una misma gota
        % TODO generalizar para varias gotas
        lmin = zeros(numnodes,1);
            for k = 1:numnodes
               % extraiga los nodos vecinos a un nodo en la misma gota 
               nodesadj = geom.nodecon2node{k};
               lmin(k) = min(normesp(repmat(geom.nodes(k,:),[size(nodesadj,1) 1]) - geom.nodes(nodesadj,:)));  
            end
    elseif geom.numdrops > 1
        % lmin entre nodos de una misma gota y entre varias gotas
        % lmin entre nodos de una misma gota
        % TODO generalizar para varias gotas
        lmin = zeros(numnodes,1);
        for j = 1:geom.numdrops
            for k = geom.nnodesdrop(j,1):geom.nnodesdrop(j,2)
               % extraiga los nodos vecinos a un nodo en la misma gota 
               nodesadj = geom.nodecon2node{k};
               lmintemp1 = min(normesp(repmat(geom.nodes(k,:),[size(nodesadj,1) 1]) - geom.nodes(nodesadj,:)));  
               % calcule para el punto respeto de los nodos de ls otras gotas
               if j == 1
               cantnodes = length(geom.nnodesdrop(2,1):geom.nnodesdrop(2,2));
               lmintemp2 = min(normesp(repmat(geom.nodes(k,:),[cantnodes 1]) - ...
                   geom.nodes(geom.nnodesdrop(2,1):geom.nnodesdrop(2,2),:)));  
               elseif j == 2
               cantnodes = length(geom.nnodesdrop(1,1):geom.nnodesdrop(1,2));
               lmintemp2 = min(normesp(repmat(geom.nodes(k,:),[cantnodes 1]) - ...
                   geom.nodes(geom.nnodesdrop(1,1):geom.nnodesdrop(1,2),:)));  
               end
               lmin(k) = min([lmintemp1 lmintemp2]);
            end
        end
    end
        % Longitud minima para verificacion de paso de tiempo.
    lmint = min(lmin);
    deltat = lmint^1.5/redfactor;
    adaptparms.lmin = lmin;
% primer paso de runge kutta
   % invoque el problema de flujo de stokes
   tic
   [velnode0,geom] = stokesmdrop(geom,parms);
   % inveoque la adaptacion de la malla
%    veladapt0 = meshadapt(geom,adaptparms,velnode0);
   
   if velopt == 1
       % hidrodinamica + adaptacion
       veladapt0 = 0;
       veladapt0 = meshadapt2(geom,adaptparms);
       k1 = deltat.*(velnode0 + veladapt0);
   elseif velopt == 2
       % normal + adaptacion
       veladapt0 = 0;
       velnormal0 = repmat(sum(velnode0.*geom.normal,2),[1 3]).*geom.normal;
       k1 = deltat.*(velnormal0 + veladapt0);
   end
   
   nodesori = geom.nodes;
   geom.nodes = geom.nodes + (1/2).*k1;
   
% segundo paso de runge kutta
    % invoque el problema de flujo de stokes
   [velnode,geom] = stokesmdrop(geom,parms);
   toc
   % inveoque la adaptacion de la malla
%    veladapt = meshadapt(geom,adaptparms,velnode);
   
   if velopt == 1
       % hidrodinamica + adaptacion
       veladapt = 0;
       veladapt = meshadapt2(geom,adaptparms);
       k2 = deltat.*(velnode + veladapt);
   elseif velopt == 2
       % normal + adaptacion
       veladapt = 0;
       velnormal = repmat(sum(velnode.*geom.normal,2),[1 3]).*geom.normal;
       k2 = deltat.*(velnormal + veladapt);
   end
  
   geom.nodes = nodesori + k2;

% escalaje
    geomprop = normalandgeo(geom,normalandgeoopt);
    geom.normalele = geomprop.normalele;
    geom.normal = geomprop.normal;
    geom.dsi = geomprop.dsi;
    geom.ds = geomprop.ds;
    geom.s = geomprop.s;
    geom.vol = geomprop.vol;
    geom.jacmat = geomprop.jacmat;
    
    errorvol = abs((geom.vol - geom.volini)./geom.volini)

    for r=1:geom.numdrops
       if errorvol(r) > errorvoltol
           % invoque escalaje
           geom = escaling(geom,optesc,r);
       end
    end

% grafique la geometria
    figure(1);
    grafscfld(geom,geom.curv);
    axis equal; view(90,0); xlabel('x1'); ylabel('x2'); zlabel('x3')
    getframe;
% error de volumen
    errorvol = abs(geom.volini - geom.vol)/geom.volini
% velocidad normal maxima y tiempo de simulacion
    velcont = max(abs(sum(velnode.*geom.normal,2)));
    geom.tiempo = geom.tiempo + deltat;
    geom.deltat = deltat;
    
% guarde resultados
    if counter == outputfreq
        itsaved = itsaved + 1;
        counter = 0;
        nombrearchivo = ['ca' num2str(ca) 'it' num2str(itsaved), '.mat'];
        save(nombrearchivo,'geom','velnode','parms','adim');
    end

end
