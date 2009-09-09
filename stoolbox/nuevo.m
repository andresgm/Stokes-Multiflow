clear;clc
% Algoritmo de flujo de stokes con sulfactantes.
ca = 0.25;
lamda = 1;
% Adimensionalizacion
adim = 1;
% Banderas de fuerza dif 0: si. 1: no
ka = 1;
kb = 0;
kc = 0;
% sulfactantes si kc = 1. 
% maranmodel = 1(lineal) definir beta y pe
% maranmodel = 2(logaritmico) definir x, e y pe
    maranmodel = 1;
    e = 0.2;
    % constante Kbar del modelo
    x = 0.75;
    % constante para modelo lineal
    beta = 0.1;
    % Peclet para la evolucion de surfactantes
    pe = 10;
% malla a usar
reflevel = 3;
% Posicion inicial de la esfera en radios
PosRadii = 0;

% pasos de tiempo de la simulacion
numtimesteps = 20000;
deltat = 0.001;
redfactor = 10;

% parametros de adaptacion
% velopt: 1 hidrodinamica velopt:2 normal
velopt = 2;
% adaptparms.opt: 4 bazhlekov
adaptparms.psi = 1;
adaptparms.lamda = lamda;
adaptparms.opt  = 4; 
%% 

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

    if kc == 1
        parms.maran.rkmaran = 1/ca;
        if maranmodel == 1
            parms.maran.maranmodel = 'linear';
            parms.maran.beta = beta;
            parms.maran.pe = pe;
        elseif maranmodel  == 2
            parms.maran.maranmodel = 'log';
            parms.maran.pe = pe;
            parms.maran.x = x;
            parms.maran.e = e;
        end
    else
        parms.maran.rkmaran = 0; 
    end
    
    
elseif adim == 2
    % adimensionalizacion de pozrkidis
    
end

%% procesamiento de la malla
% cargue el archivo base
filename = ['sph ref ' num2str(reflevel) '.mat'];
load([cd '/' filename]);
% PROCESAMIENTO DE LA MALLA ORIGINAL
% Volumen original
Radius = max(normesp(Nodes));
% Numero de Elementos y numero de Nodos
geom.nodes = Nodes;
geom.elements = Elements;
geom.numnodes = size(geom.nodes,1);
% Posicione la esfera a su punto de partida en radios 
geom.nodes(:,3) = geom.nodes(:,3) + ...
    repmat(PosRadii*Radius,geom.numnodes,1);
geom.numelements = size(Elements,1);
% Index table
geom.indextable = [1:1:geom.numnodes];
    % Tabla de elementos singulares a cada nodo
geom.element2node = element2node(geom.elements);
    % Elementos singulares a cada nodo tabla SPARSE
geom.e2nsparse = ele2nodesp(geom.elements);
    % Tabla de conectividad de nodos, bordes, e.t.c
[geom.nodecon2node,geom.edgeindex] = ...
    node2node(geom.elements,geom.element2node,geom.indextable);   
    
% calcule el volumen inicial de la gota
normalandgeoopt.normal = 1;
normalandgeoopt.areas = 1;
normalandgeoopt.vol = 1;
geom = normalandgeo(geom,normalandgeoopt)
geom.volini = geom.vol;

%% evolucion mediante runge kutta de 4to orden
counter = 0;
for p=1:numtimesteps
% calcule la distancia minima de adaptacion y el paso de tiempo
    counter = counter + 1;    
    l = zeros(size(geom.edgeindex,1),1);
    for k=1:size(geom.edgeindex,1)
        % Extraiga los nodos del kEdge
        kedge = geom.edgeindex(k,:);
        % Longitud del kEdge
        l(k) = normesp(geom.nodes(kedge(1),:) - geom.nodes(kedge(2),:));
    end
    % Longitud minima para verificacion de paso de tiempo.
    lmin = min(l);
    adaptparms.lmin = lmin;

    deltat = lmin^1.5/redfactor;
    deltat = 0.01;
% primer paso de runge kutta
   % invoque el problema de flujo de stokes
   [velnode0,geom] = stokes(geom,parms);
   % inveoque la adaptacion de la malla
   veladapt0 = meshadapt(geom,adaptparms,velnode0);
   
   if velopt == 1
       % hidrodinamica + adaptacion
       k1 = deltat.*(velnode0 + veladapt0);
   elseif velopt == 2
       % normal + adaptacion
       velnormal0 = repmat(sum(velnode0.*geom.normal,2),[1 3]).*geom.normal;
       k1 = deltat.*(velnormal0 + veladapt0);
   end
   
   nodeori = geom.nodes;
   geom.nodes = geom.nodes + (1/2).*k1;
   
% segundo paso de runge kutta
    % invoque el problema de flujo de stokes
   [velnode,geom] = stokes(geom,parms);
   % inveoque la adaptacion de la malla
   veladapt = meshadapt(geom,adaptparms,velnode);
   
   if velopt == 1
       % hidrodinamica + adaptacion
       k2 = deltat.*(velnode + veladapt);
   elseif velopt == 2
       % normal + adaptacion
       velnormal = repmat(sum(velnode.*geom.normal,2),[1 3]).*geom.normal;
       k2 = deltat.*(velnormal + veladapt);
   end
   
   nodesori = geom.nodes;
   
   geom.nodes = nodesori + k2;
   
% grafique la geometria
 grafscfld(geom,geom.curv);
 axis equal; view(90,0); xlabel('x1'); ylabel('x2'); zlabel('x3')
 getframe;
% error de volumen
errorvol = abs(geom.volini - geom.vol)/geom.volini
% velocidad normal maxima
velcont = max(abs(sum(velnode.*geom.normal,2)))
end
