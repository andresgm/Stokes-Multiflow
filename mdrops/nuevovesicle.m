% CALCULO DEL FLUJO DE STOKES PARA UNA VESICULA
% IMPLEMENTADO FLUJO INFINITO, SEMIINFINITO
% IMPLEMENTADO SINGLE Y DOUBLE LAYER
clear;clc;%close all;
%% opciones de carga de archivos
    % nombre de archivo a cargar y carpeta
nombreorigen = 'sph ref 3';
carpetaorigen = '';
iteracion = [];

    % nombre de archivo a guardar y carpeta
nombredestino = 'it';
carpetadestino = 'sedimentacion_vesicle_g0_100_kbar20_electrostatic';
    % simulacion nueva desde cero optsim = 0
    % continue la simulacion optsim = 1
    % simulacion nueva desde archivo de resultados optsim = 2
opcionsim = 0;

% Algoritmo de flujo de stokes con surfactantes.
ca = 0;
lamda = 1;
g0 = 100;
e0 = 0;
% tipo de flujo flow: 'inf'  flow:'semiinf'
flow = 'semiinf';
% opcion de calculo de la curvatura 1: paraboloid fitting; 2: extended par;
% 3: basado en laplace beltrami
curvopt = 3;
% Constantes del modelo de bending
    % Constante c del modelo
c = 0.1;
    % Coeficiente de rigides al doblamiento en KbT (kappa)
kbar = 20;
    % Coeficiente adimensional de resistencia al cambio de area: Ka*R_0^2/kappa.
    % Usar un numero entre 2e2 y 9e5? %TODO
kext = 2.5e4;


% Adimensionalizacion
adim = 1;
% frecuencia de guardar resultados
outputfreq = 10;

% Banderas de fuerza dif 0: si. 1: no
    % Tension (superficial)
ka = 1;
    % gravedad
kb = 1;
    % bending
kc = 1;
    % campo electrico
kd = 0;
    % interaccion electrostatica
ke = 1;


lie = 48.75;
psi1ie = 47.73;
psi2ie = 3.36;
gammaie = 26300;


% numero de gotas
geom.numdrops = 1;
% Coordenadas de los centroides de las gotas
xc =[0 0 10];
% Introduzca el/los radios de la/s gotas
xr=[1];

% pasos de tiempo de la simulacion
numtimesteps = 80000;

redfactor = 10;

% Tipo de integracion 1:Runge Kutta segundo orden 2:Runge Kutta cuarto orden
% 3: Adams-Bashford
inttype = 3;

% parametros de adaptacion
% velopt: 1 hidrodinamica velopt:2 normal velopt:3 passive (zinchenko et al.)
velopt = 3;
% meshadapt lowenberg
adaptparms.psi = 1;
adaptparms.lamda = lamda;
% escalaje
errorvoltol = 1e-6;
optesc.maxit = 15000;
optesc.kp = 20;
optesc.deltate = 0.01;
optesc.tolerrorvol = errorvoltol;

%% procesamiento de parametros
if adim == 1
    
    warning...
 ('DeltaF = (Sigma/g0)(2H) + Z - (1/g0)(4H^3 + 2Laps(H) + 4KH)');
    
    parms.flow = flow;
    parms.w = 0;
    % adimensionalizacion del single layer
    parms.rkextf = 2*ca/g0;
    parms.rksl = 2;
    parms.rkdl = 2*(lamda - 1)/(lamda + 1);

    parms.g0 = g0;   
    parms.lamda = lamda;
    parms.ca = ca;
    
    % parametros de simulacion
    if ka == 1 
        % curvatura constante
        parms.rkcurv = 1/g0;
    else
        parms.rkcurv = 0;
    end

    if kb == 1
        % gravedad
        parms.rkgrav = 1;
    else
        parms.rkgrav = 0;
    end    
    
    if kc == 1
        % bending
        parms.rkbend = 1/g0;
        parms.bending.c = c;
        parms.bending.kbar = kbar;
        parms.bending.kext = kext;
        parms.bending.sigma = 0;
    else
        parms.rkbend = 0;
    end    
    
    if kd == 1
        % campo electrico
        parms.rkelect = e0/g0;
    else
        parms.rkelect = 0;
    end
    
    if ke == 1
       % Interaccion electrostatica
       parms.rkelestat = gammaie;
       parms.elestat.l = lie;
       parms.elestat.psi1 = psi1ie;
       parms.elestat.psi2 = psi2ie;
    end
    
else
    error('Wrong adim parameter')
end
parms.curvopt = curvopt;
% guarde temporalmente los parametros
parmstemp = parms;
%% procesamiento de la malla
sbar = systembar();

if opcionsim == 0
    % cargue el archivo base
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
    geom = drops(geom,xc,xr);
    
    if parms.lamda ~= 1
       geom.W = zeros(numnodes,3); 
       geom.velnodeant = zeros(numnodes,3);
    end

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
    
    % Encuentre los vertices de la malla si se va a usar la adaptacion de malla
    % pasiva de Zinchenco et al. 1997
    if velopt == 3
       geom.vertices = extractvertices(geom);
    end

    % calcule el volumen inicial de la gota
    normalandgeoopt.normal = 1;
    normalandgeoopt.areas = 1;
    normalandgeoopt.vol = 1;
    geomprop = normalandgeo(geom,normalandgeoopt);
%     geom.normalele = geomprop.normalele;
    geom.normal = geomprop.normal;
    geom.dsi = geomprop.dsi;
    geom.ds = geomprop.ds;
    geom.s = geomprop.s;
    geom.vol = geomprop.vol;
    geom.jacmat = geomprop.jacmat;
    geom.volini = geom.vol;
    geom.areaini = geom.s;
    geom.xcini = centroide(geom);
    
    % Calculo inicial de la curvatura usando ajuste a superficie cuadratica
    
    paropt.tipo = 'extended';
    [geom.curv,geom.normal,geom.Kg] = curvparaboloid(geom,paropt);
    
    if isempty(carpetadestino) == 1
        direccion = [cd  sbar nombredestino num2str(iteracion) '.mat'];   
    else
        direccion = ...
         [cd  sbar carpetadestino sbar nombredestino num2str(iteracion) '.mat'];        
    end
    
    direcciondestino = [cd sbar carpetadestino sbar nombredestino];
    mkdir([cd sbar,carpetadestino]);
        
    paso = 1;
    counter = 0;
    geom.tiempo = 0;
    itsaved = 0;

elseif opcionsim == 1
    % cargue desde resultados y continue la simulacion
    carpetadestino = carpetaorigen;
    nombredestino = nombreorigen;
    if isempty(carpetaorigen) == 1
        direccion = [cd  sbar nombreorigen num2str(iteracion) '.mat'];        
    else
        direccion = ...
         [cd  sbar carpetaorigen sbar nombreorigen num2str(iteracion) '.mat'];
    end
    
    load(direccion);
    
    paso = iteracion + 1;
    counter = 0;
    itsaved = iteracion;

    direcciondestino = [cd  sbar carpetadestino sbar nombredestino];
    mkdir([cd sbar,carpetadestino]);
    numnodes = size(geom.nodes,1);
    numelements = size(geom.elements,1);    
    normalandgeoopt.normal = 1;
    normalandgeoopt.areas = 1;
    normalandgeoopt.vol = 1;
    
    parms.curvopt = parmstemp.curvopt;
elseif opcionsim == 2
    % cargue desde resultados y realice una nueva simulacion
    % cargue desde resultados y continue la simulacion
    if isempty(carpetaorigen) == 1
        direccion = [cd  sbar nombreorigen num2str(iteracion) '.mat'];        
        
    else
        direccion = ...
         [cd  sbar carpetaorigen sbar nombreorigen num2str(iteracion) '.mat'];        
    end
    
    load(direccion);
    
    direcciondestino = [cd  sbar carpetadestino sbar nombredestino];
    mkdir([cd sbar,carpetadestino]);
        
    paso = 1;    
    if parmstemp.rkbend ~= 0
        % asigne el bending de la simulacion cargada
        parmstemp.bending.sigma = parms.bending.sigma;
    end
    parms = parmstemp;
    counter = 0;
    geom.tiempo = 0;
    itsaved = 0;
    numnodes = size(geom.nodes,1);
    numelements = size(geom.elements,1);
    normalandgeoopt.normal = 1;
    normalandgeoopt.areas = 1;
    normalandgeoopt.vol = 1;
end


%% Ciclo principal
xcant = centroide(geom);
geom.velcentroid = [0 0 0];
geom.xc = xcant;
if velopt == 3
   veladapt = zeros(numnodes,3);
end
for p = paso:numtimesteps
tic
% calcule la distancia minima de adaptacion y el paso de tiempo
    counter = counter + 1;
    disp(['iteracion = ', num2str(p)])
    
%     if p == 4
%        profile on
%     end
   
    if geom.numdrops == 1
        % lmin entre nodos de una misma gota
        % TODO generalizar para varias gotas
        lmin = zeros(numnodes,1);
            for k = 1:numnodes
               % extraiga los nodos vecinos a un nodo en la misma gota 
               nodesadj = geom.nodecon2node{k};
               lmin(k) = min(normesp(repmat(geom.nodes(k,:),[size(nodesadj,1) 1])...
                  - geom.nodes(nodesadj,:)));  
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
               lmintemp1 = min(normesp(repmat(geom.nodes(k,:),[size(nodesadj,1) 1])...
                  - geom.nodes(nodesadj,:)));  
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
    parms.lmin = lmin;
    
    if inttype == 1
       % Runge-Kutta de segundo orden
       
       % primer paso de runge kutta f1
       % invoque el problema de flujo de stokes
       [velnode1,geom,parms] = stokesvesicle(geom,parms);
       % invoque la adaptacion de la malla
       % veladapt0 = meshadapt(geom,adaptparms,velnode0);

       if velopt == 1
           % hidrodinamica + adaptacion
           veladapt1 = meshadapt2(geom,parms);
           temporal = repmat(geom.velcentroid,[numnodes 1]);
    %        veltan = temporal - repmat(sum(temporal.*geom.normal,2),[1 3]).*geom.normal;
           veltan = 0;
           f1 = (velnode1 + veladapt1 + veltan);
       elseif velopt == 2
           % normal + adaptacion
           temporal = repmat(geom.velcentroid,[numnodes 1]);
           veltan = temporal - repmat(sum(temporal.*geom.normal,2),[1 3]).*geom.normal;
    %        veltan = 0;
    %        veladapt = 0;
           velnormal1 = repmat(sum(velnode1.*geom.normal,2),[1 3]).*geom.normal;
           veladapt1 = meshadapt2(geom,parms);
           f1 = (velnormal1 + veladapt1 + veltan);
       elseif velopt == 3
           % passive (zinchenco et al. 1997)
           velnormal = repmat(sum(velnode1.*geom.normal,2),[1 3]).*geom.normal;
           veladapt = meshadaptgrad(geom,velnormal,veladapt);
           f1 = (velnormal + veladapt);
       end

       nodes0 = geom.nodes;
       geom.nodes = geom.nodes + (1/2)*deltat*f1;
       
       paropt.tipo = 'extended';
       [dummyc,geom.normal,dummyk] = curvparaboloid(geom,paropt);
       
       % segundo paso de runge kutta f4
       % invoque el problema de flujo de stokes
       [velnode,geom,parms] = stokesvesicle(geom,parms);
       % invoque la adaptacion de la malla
    %    veladapt0 = meshadapt(geom,adaptparms,velnode0);

       if velopt == 1
           % hidrodinamica + adaptacion
           veladapt = meshadapt2(geom,parms);
           temporal = repmat(geom.velcentroid,[numnodes 1]);
    %        veltan = temporal - repmat(sum(temporal.*geom.normal,2),[1 3]).*geom.normal;
           veltan = 0;
           f2 = (velnode + veladapt + veltan);
       elseif velopt == 2
           % normal + adaptacion
           temporal = repmat(geom.velcentroid,[numnodes 1]);
           veltan = temporal - repmat(sum(temporal.*geom.normal,2),[1 3]).*geom.normal;
    %        veltan = 0;
    %        veladapt = 0;
           velnormal = repmat(sum(velnode.*geom.normal,2),[1 3]).*geom.normal;
           veladapt = meshadapt2(geom,parms);
           f2 = (velnormal + veladapt + veltan);
       elseif velopt == 3
           % passive (zinchenco et al. 1997)
           velnormal = repmat(sum(velnode.*geom.normal,2),[1 3]).*geom.normal;
           veladapt = meshadaptgrad(geom,velnormal,veladapt);
           f2 = (velnormal + veladapt);
       end
       
       geom.nodes = nodes0 + deltat*f2;
       
       paropt.tipo = 'extended';
       [dummyc,geom.normal,dummyk] = curvparaboloid(geom,paropt);
       
    elseif inttype == 2
       % Runge-Kutta de cuarto orden
       
       % primer paso de runge kutta f1
       % invoque el problema de flujo de stokes
       [velnode1,geom,parms] = stokesvesicle(geom,parms);
       % invoque la adaptacion de la malla
       % veladapt0 = meshadapt(geom,adaptparms,velnode0);

       if velopt == 1
           % hidrodinamica + adaptacion
           veladapt1 = meshadapt2(geom,parms);
           temporal = repmat(geom.velcentroid,[numnodes 1]);
    %        veltan = temporal - repmat(sum(temporal.*geom.normal,2),[1 3]).*geom.normal;
           veltan = 0;
           f1 = (velnode1 + veladapt1 + veltan);
       elseif velopt == 2
           % normal + adaptacion
           temporal = repmat(geom.velcentroid,[numnodes 1]);
           veltan = temporal - repmat(sum(temporal.*geom.normal,2),[1 3]).*geom.normal;
    %        veltan = 0;
    %        veladapt = 0;
           velnormal1 = repmat(sum(velnode1.*geom.normal,2),[1 3]).*geom.normal;
           veladapt1 = meshadapt2(geom,parms);
           f1 = (velnormal1 + veladapt1 + veltan);
       elseif velopt == 3
           % passive (zinchenco et al. 1997)
           velnormal = repmat(sum(velnode1.*geom.normal,2),[1 3]).*geom.normal;
           veladapt = meshadaptgrad(geom,velnormal,veladapt);
           f1 = (velnormal + veladapt);
       end

       nodes0 = geom.nodes;

    %% segundo paso de runge kutta f2
       % invoque el problema de flujo de stokes
       geom.nodes = geom.nodes + (1/2)*deltat*f1;
       
       paropt.tipo = 'extended';
       [dummyc,geom.normal,dummyk] = curvparaboloid(geom,paropt);
       
       [velnode2,geom,parms] = stokesvesicle(geom,parms);
       % invoque la adaptacion de la malla
    %    veladapt0 = meshadapt(geom,adaptparms,velnode0);

       if velopt == 1
           % hidrodinamica + adaptacion
           veladapt2 = meshadapt2(geom,parms);
           temporal = repmat(geom.velcentroid,[numnodes 1]);
    %        veltan = temporal - repmat(sum(temporal.*geom.normal,2),[1 3]).*geom.normal;
           veltan = 0;
           f2 = (velnode2 + veladapt2 + veltan);
       elseif velopt == 2
           % normal + adaptacion
           temporal = repmat(geom.velcentroid,[numnodes 1]);
           veltan = temporal - repmat(sum(temporal.*geom.normal,2),[1 3]).*geom.normal;
    %        veltan = 0;
    %        veladapt = 0;
           velnormal2 = repmat(sum(velnode2.*geom.normal,2),[1 3]).*geom.normal;
           veladapt2 = meshadapt2(geom,parms);
           f2 = (velnormal2 + veladapt2 + veltan);
       elseif velopt == 3
           % passive (zinchenco et al. 1997)
           velnormal = repmat(sum(velnode2.*geom.normal,2),[1 3]).*geom.normal;
           veladapt = meshadaptgrad(geom,velnormal,veladapt);
           f2 = (velnormal + veladapt);
       end
       
    %% tercer paso de runge kutta f3
       % invoque el problema de flujo de stokes
       geom.nodes = nodes0 + (1/2)*deltat*f2;
       
       paropt.tipo = 'extended';
       [dummyc,geom.normal,dummyk] = curvparaboloid(geom,paropt);
       
       [velnode3,geom,parms] = stokesvesicle(geom,parms);
       % invoque la adaptacion de la malla
    %    veladapt0 = meshadapt(geom,adaptparms,velnode0);

       if velopt == 1
           % hidrodinamica + adaptacion
           veladapt3 = meshadapt2(geom,parms);
           temporal = repmat(geom.velcentroid,[numnodes 1]);
    %        veltan = temporal - repmat(sum(temporal.*geom.normal,2),[1 3]).*geom.normal;
           veltan = 0;
           f3 = (velnode3 + veladapt3 + veltan);
       elseif velopt == 2
           % normal + adaptacion
           temporal = repmat(geom.velcentroid,[numnodes 1]);
           veltan = temporal - repmat(sum(temporal.*geom.normal,2),[1 3]).*geom.normal;
    %        veltan = 0;
    %        veladapt = 0;
           velnormal3 = repmat(sum(velnode3.*geom.normal,2),[1 3]).*geom.normal;
           veladapt3 = meshadapt2(geom,parms);
           f3 = (velnormal3 + veladapt3 + veltan);
       elseif velopt == 3
           % passive (zinchenco et al. 1997)
           velnormal = repmat(sum(velnode3.*geom.normal,2),[1 3]).*geom.normal;
           veladapt = meshadaptgrad(geom,velnormal,veladapt);
           f3 = (velnormal + veladapt);
       end
       
    %% cuarto paso de runge kutta f4
       % invoque el problema de flujo de stokes
       geom.nodes = nodes0 + deltat*f3;
       
       paropt.tipo = 'extended';
       [dummyc,geom.normal,dummyk] = curvparaboloid(geom,paropt);
       
       [velnode,geom,parms] = stokesvesicle(geom,parms);
       % invoque la adaptacion de la malla
    %    veladapt0 = meshadapt(geom,adaptparms,velnode0);

       if velopt == 1
           % hidrodinamica + adaptacion
           veladapt = meshadapt2(geom,parms);
           temporal = repmat(geom.velcentroid,[numnodes 1]);
    %        veltan = temporal - repmat(sum(temporal.*geom.normal,2),[1 3]).*geom.normal;
           veltan = 0;
           f4 = (velnode + veladapt + veltan);
       elseif velopt == 2
           % normal + adaptacion
           temporal = repmat(geom.velcentroid,[numnodes 1]);
           veltan = temporal - repmat(sum(temporal.*geom.normal,2),[1 3]).*geom.normal;
    %        veltan = 0;
    %        veladapt = 0;
           velnormal = repmat(sum(velnode.*geom.normal,2),[1 3]).*geom.normal;
           veladapt = meshadapt2(geom,parms);
           f4 = (velnormal + veladapt + veltan);
       elseif velopt == 3
           % passive (zinchenco et al. 1997)
           velnormal = repmat(sum(velnode.*geom.normal,2),[1 3]).*geom.normal;
           veladapt = meshadaptgrad(geom,velnormal,veladapt);
           f4 = (velnormal + veladapt);
       end

       geom.nodes = nodes0 + deltat*(f1+2*f2+2*f3+f4)/6;
       
       paropt.tipo = 'extended';
       [dummyc,geom.normal,dummyk] = curvparaboloid(geom,paropt);
       
   elseif inttype == 3
    
   %% Para los primeros pasos usamos RK4 para inicializar los puntos de
   %% Adams-Bashforth-Moulton
      if p <=3
          % primer paso de runge kutta f1
          % invoque el problema de flujo de stokes
          [velnode1,geom,parms] = stokesvesicle(geom,parms);
          % invoque la adaptacion de la malla
          % veladapt0 = meshadapt(geom,adaptparms,velnode0);

          if velopt == 1
              % hidrodinamica + adaptacion
              veladapt1 = meshadapt2(geom,parms);
              temporal = repmat(geom.velcentroid,[numnodes 1]);
       %        veltan = temporal - repmat(sum(temporal.*geom.normal,2),[1 3]).*geom.normal;
              veltan = 0;
              f1 = (velnode1 + veladapt1 + veltan);
          elseif velopt == 2
              % normal + adaptacion
              temporal = repmat(geom.velcentroid,[numnodes 1]);
              veltan = temporal - repmat(sum(temporal.*geom.normal,2),[1 3]).*geom.normal;
       %        veltan = 0;
       %        veladapt = 0;
              velnormal1 = repmat(sum(velnode1.*geom.normal,2),[1 3]).*geom.normal;
              veladapt1 = meshadapt2(geom,parms);
              f1 = (velnormal1 + veladapt1 + veltan);
          elseif velopt == 3
              % passive (zinchenco et al. 1997)
              velnormal = repmat(sum(velnode1.*geom.normal,2),[1 3]).*geom.normal;
              veladapt = meshadaptgrad(geom,velnormal,veladapt);
              f1 = (velnormal + veladapt);
          end

          abm{p} = f1;
          nodes0 = geom.nodes;

       %% segundo paso de runge kutta f2
          % invoque el problema de flujo de stokes
          geom.nodes = geom.nodes + (1/2)*deltat*f1;
          
          paropt.tipo = 'extended';
          [dummyc,geom.normal,dummyk] = curvparaboloid(geom,paropt);
          
          [velnode2,geom,parms] = stokesvesicle(geom,parms);
          % invoque la adaptacion de la malla
       %    veladapt0 = meshadapt(geom,adaptparms,velnode0);

          if velopt == 1
              % hidrodinamica + adaptacion
              veladapt2 = meshadapt2(geom,parms);
              temporal = repmat(geom.velcentroid,[numnodes 1]);
       %        veltan = temporal - repmat(sum(temporal.*geom.normal,2),[1 3]).*geom.normal;
              veltan = 0;
              f2 = (velnode2 + veladapt2 + veltan);
          elseif velopt == 2
              % normal + adaptacion
              temporal = repmat(geom.velcentroid,[numnodes 1]);
              veltan = temporal - repmat(sum(temporal.*geom.normal,2),[1 3]).*geom.normal;
       %        veltan = 0;
       %        veladapt = 0;
              velnormal2 = repmat(sum(velnode2.*geom.normal,2),[1 3]).*geom.normal;
              veladapt2 = meshadapt2(geom,parms);
              f2 = (velnormal2 + veladapt2 + veltan);
          elseif velopt == 3
              % passive (zinchenco et al. 1997)
              velnormal = repmat(sum(velnode2.*geom.normal,2),[1 3]).*geom.normal;
              veladapt = meshadaptgrad(geom,velnormal,veladapt);
              f2 = (velnormal + veladapt);
          end

       %% tercer paso de runge kutta f3
          % invoque el problema de flujo de stokes
          geom.nodes = nodes0 + (1/2)*deltat*f2;
          
          paropt.tipo = 'extended';
          [dummyc,geom.normal,dummyk] = curvparaboloid(geom,paropt);
          
          [velnode3,geom,parms] = stokesvesicle(geom,parms);
          % invoque la adaptacion de la malla
       %    veladapt0 = meshadapt(geom,adaptparms,velnode0);

          if velopt == 1
              % hidrodinamica + adaptacion
              veladapt3 = meshadapt2(geom,parms);
              temporal = repmat(geom.velcentroid,[numnodes 1]);
       %        veltan = temporal - repmat(sum(temporal.*geom.normal,2),[1 3]).*geom.normal;
              veltan = 0;
              f3 = (velnode3 + veladapt3 + veltan);
          elseif velopt == 2
              % normal + adaptacion
              temporal = repmat(geom.velcentroid,[numnodes 1]);
              veltan = temporal - repmat(sum(temporal.*geom.normal,2),[1 3]).*geom.normal;
       %        veltan = 0;
       %        veladapt = 0;
              velnormal3 = repmat(sum(velnode3.*geom.normal,2),[1 3]).*geom.normal;
              veladapt3 = meshadapt2(geom,parms);
              f3 = (velnormal3 + veladapt3 + veltan);
          elseif velopt == 3
              % passive (zinchenco et al. 1997)
              velnormal = repmat(sum(velnode3.*geom.normal,2),[1 3]).*geom.normal;
              veladapt = meshadaptgrad(geom,velnormal,veladapt);
              f3 = (velnormal + veladapt);
          end

       %% cuarto paso de runge kutta f4
          % invoque el problema de flujo de stokes
          geom.nodes = nodes0 + deltat*f3;
          
          paropt.tipo = 'extended';
          [dummyc,geom.normal,dummyk] = curvparaboloid(geom,paropt);
          
          [velnode,geom,parms] = stokesvesicle(geom,parms);
          % invoque la adaptacion de la malla
       %    veladapt0 = meshadapt(geom,adaptparms,velnode0);

          if velopt == 1
              % hidrodinamica + adaptacion
              veladapt = meshadapt2(geom,parms);
              temporal = repmat(geom.velcentroid,[numnodes 1]);
       %        veltan = temporal - repmat(sum(temporal.*geom.normal,2),[1 3]).*geom.normal;
              veltan = 0;
              f4 = (velnode + veladapt + veltan);
          elseif velopt == 2
              % normal + adaptacion
              temporal = repmat(geom.velcentroid,[numnodes 1]);
              veltan = temporal - repmat(sum(temporal.*geom.normal,2),[1 3]).*geom.normal;
       %        veltan = 0;
       %        veladapt = 0;
              velnormal = repmat(sum(velnode.*geom.normal,2),[1 3]).*geom.normal;
              veladapt = meshadapt2(geom,parms);
              f4 = (velnormal + veladapt + veltan);
          elseif velopt == 3
              % passive (zinchenco et al. 1997)
              velnormal = repmat(sum(velnode.*geom.normal,2),[1 3]).*geom.normal;
              veladapt = meshadaptgrad(geom,velnormal,veladapt);
              f4 = (velnormal + veladapt);
          end

          geom.nodes = nodes0 + deltat*(f1+2*f2+2*f3+f4)/6;
          
          paropt.tipo = 'extended';
          [dummyc,geom.normal,dummyk] = curvparaboloid(geom,paropt);

      else

          % Ahora podemos usar el m?todo Predictor-Corrector de
          % Adams-Bashforth-Moulton

          % Paso predictor

          % Calculo de la velocidad en el punto actual

          [velnode1,geom,parms] = stokesvesicle(geom,parms);
          % invoque la adaptacion de la malla
          % veladapt0 = meshadapt(geom,adaptparms,velnode0);

          if velopt == 1
              % hidrodinamica + adaptacion
              veladapt1 = meshadapt2(geom,parms);
              temporal = repmat(geom.velcentroid,[numnodes 1]);
       %        veltan = temporal - repmat(sum(temporal.*geom.normal,2),[1 3]).*geom.normal;
              veltan = 0;
              f1 = (velnode1 + veladapt1 + veltan);
          elseif velopt == 2
              % normal + adaptacion
              temporal = repmat(geom.velcentroid,[numnodes 1]);
              veltan = temporal - repmat(sum(temporal.*geom.normal,2),[1 3]).*geom.normal;
       %        veltan = 0;
       %        veladapt = 0;
              velnormal1 = repmat(sum(velnode1.*geom.normal,2),[1 3]).*geom.normal;
              veladapt1 = meshadapt2(geom,parms);
              f1 = (velnormal1 + veladapt1 + veltan);
          elseif velopt == 3
              % passive (zinchenco et al. 1997)
              velnormal = repmat(sum(velnode1.*geom.normal,2),[1 3]).*geom.normal;
              veladapt = meshadaptgrad(geom,velnormal,veladapt);
              f1 = (velnormal + veladapt);
          end

          abm{4} = f1;
          nodes0 = geom.nodes;
          geom.nodes = geom.nodes + deltat*(-9*abm{1}+37*abm{2}-59*abm{3}+55*abm{4})/24;
          
          paropt.tipo = 'extended';
          [dummyc,geom.normal,dummyk] = curvparaboloid(geom,paropt);

          % Paso corrector

          % Calculo de la velocidad en el punto siguiente con la prediccion

          [velnode,geom,parms] = stokesvesicle(geom,parms);
          % invoque la adaptacion de la malla
          % veladapt0 = meshadapt(geom,adaptparms,velnode0);

          if velopt == 1
              % hidrodinamica + adaptacion
              veladapt = meshadapt2(geom,parms);
              temporal = repmat(geom.velcentroid,[numnodes 1]);
       %        veltan = temporal - repmat(sum(temporal.*geom.normal,2),[1 3]).*geom.normal;
              veltan = 0;
              f2 = (velnode + veladapt + veltan);
          elseif velopt == 2
              % normal + adaptacion
              temporal = repmat(geom.velcentroid,[numnodes 1]);
              veltan = temporal - repmat(sum(temporal.*geom.normal,2),[1 3]).*geom.normal;
       %        veltan = 0;
       %        veladapt = 0;
              velnormal = repmat(sum(velnode.*geom.normal,2),[1 3]).*geom.normal;
              veladapt = meshadapt2(geom,parms);
              f2 = (velnormal + veladapt + veltan);
          elseif velopt == 3
              % passive (zinchenco et al. 1997)
              velnormal = repmat(sum(velnode.*geom.normal,2),[1 3]).*geom.normal;
              veladapt = meshadaptgrad(geom,velnormal,veladapt);
              f2 = (velnormal + veladapt);
          end

          abm{5} = f2;

          geom.nodes = nodes0 + deltat*(abm{2}-5*abm{3}+19*abm{4}+9*abm{5})/24;
          
          paropt.tipo = 'extended';
          [dummyc,geom.normal,dummyk] = curvparaboloid(geom,paropt);

          % Actualizaci?n nodos

          abm{1} = abm{2};
          abm{2} = abm{3};
          abm{3} = abm{4};
          
      end
    end
%   parms.bending.sigma
%% escalaje
    geomprop = normalandgeo(geom,normalandgeoopt);
%     geom.normalele = geomprop.normalele;
    % geom.normal = geomprop.normal;
    geom.dsi = geomprop.dsi;
    geom.ds = geomprop.ds;
    geom.s = geomprop.s;
    geom.vol = geomprop.vol;
%     geom.jacmat = geomprop.jacmat;
    
    errorvol = abs((geom.vol - geom.volini)./geom.volini);

    for r=1:geom.numdrops
       if errorvol(r) > errorvoltol
           % invoque escalaje
           disp('Realizando Escalaje')
           geom = escaling(geom,optesc,r,errorvol);
       end
    end

% error de volumen
    errorvol = abs(geom.volini - geom.vol)/geom.volini;
% velocidad normal maxima y tiempo de simulacion
    velcont = max(abs(sum(velnode.*geom.normal,2)));
    geom.tiempo = geom.tiempo + deltat;
    geom.deltat = deltat;
% calculo del centroide de la gota y velocidad del centroide
    xcant = geom.xc;
    geom.xc = centroide(geom);
    geom.velcentroid = (geom.xc - xcant)./deltat;
    
% Visualizacion
    figure(1);
    grafscfld(geom,geom.curv);
    axis equal; view(90,0); xlabel('x1'); ylabel('x2'); zlabel('x3'); colorbar;
    getframe; title('curv');
    
    figure(2); plot(geom.tiempo,geom.velcentroid(:,3),'*');hold on; title('vc x3');

    figure(3); plot(geom.tiempo,geom.xc(:,3),'*');hold on; title('xc x3');    
    
    if kc ~= 0
      figure(4); grafscfld(geom,geom.lapcurv);
      axis equal; view(90,0); xlabel('x1'); ylabel('x2'); zlabel('x3'); colorbar;
      getframe; title('laplace curv');
      figure(5); plot(geom.tiempo,parms.bending.sigma,'*'); hold on;title('sigma');
    end
    
    if ke ~= 0
       figure(6); plot(geom.tiempo,max(geom.deltafelestat),'*r');hold on;
       title('Electrostat vs. Grav');
       plot(geom.tiempo,max(geom.deltafgrav),'*b')
    end
    disp(['Grav: ', num2str(max(geom.deltafgrav)),' Electroestatica: ', ...
       max(geom.deltafelestat)]);

% guarde resultados
    if counter == outputfreq
        itsaved = itsaved + 1;
        counter = 0;
        nombrearchivo = [direcciondestino num2str(itsaved), '.mat'];
        save(nombrearchivo,'geom','velnode','parms','adim','');
    end
    disp(carpetadestino)
    disp(['deltat: ', num2str(geom.deltat)])
    
%     if p == 14
%        profile viewer
%     end
    
toc
end
