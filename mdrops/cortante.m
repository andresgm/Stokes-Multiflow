% CALCULO DEL FLUJO DE STOKES PARA UNA VESICULA
% IMPLEMENTADO FLUJO INFINITO, SEMIINFINITO
% IMPLEMENTADO SINGLE Y DOUBLE LAYER
clear;clc; %close all;
%% opciones de carga de archivos
    % nombre de archivo a cargar y carpeta
nombreorigen = 'sph ref 3.mat';
carpetaorigen = '';
iteracion = [];

    % nombre de archivo a guardar y carpeta
nombredestino = 'it';
carpetadestino = 'pruebacortante_fem';
    % simulacion nueva desde cero optsim = 0
    % continue la simulacion optsim = 1
    % simulacion nueva desde archivo de resultados optsim = 2
opcionsim = 0;

% Parametros introduccion ruido para minimizar efecto de la simetria de la
% malla

noiseint = 0.05;
noiserep = 4;

% Algoritmo de flujo de stokes.
ca = 0.1;
lamda = 6.4;

% tipo de flujo flow: 'inf'  flow:'semiinf'
flow = 'inf';
% opcion de calculo de la curvatura 1: paraboloid fitting; 2: extended par;
% 3: basado en laplace beltrami
curvopt = 3;
% Constantes del modelo de bending
% Constante c del modelo
c = 0.1;
    % Coeficiente de rigides al doblamiento en KbT (kappa)
kbar = 20;
    % Coeficiente adimensional de resistencia al cambio de area:
    % Ka*R_0^2/kappa.
kext = 4e8;


% Adimensionalizacion
adim = 1;
% frecuencia de guardar resultados
outputfreq = 10;

% pasos de tiempo de la simulacion
numtimesteps = 80000;

redfactor = 10;

% Tipo de integracion 1:Runge Kutta segundo orden 2:Runge Kutta cuarto orden
% 3: Adams-Bashford
inttype = 3;

% Estamos usando solo la adaptacion de malla pasiva propuesta por Zinchenko
% et al. 1997 y 1999.
% Sin adaptacion de malla. OJO!

% escalaje
errorvoltol = 1e-6;
optesc.maxit = 15000;
optesc.kp = 20;
optesc.deltate = 0.01;
optesc.tolerrorvol = errorvoltol;

%% procesamiento de parametros
parms.flow = flow;
parms.w = 0;
% adimensionalizacion del single layer
parms.rkextf = 2*ca;
parms.rksl = 2;
parms.rkdl = 2*(lamda - 1)/(lamda + 1);
parms.lamda = lamda;
parms.ca = ca;
    
% Coeficiente termino de curvatura
parms.rkcurv = 1;

% Coeficiente adimensional termino bending
parms.rkbend = 1;
parms.bending.c = c;
parms.bending.kbar = kbar;
parms.bending.kext = kext;
parms.bending.sigma = 0; 
    
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
    
    if parms.lamda ~= 1
       geom.W = zeros(numnodes,3); 
       geom.velnodeant = zeros(numnodes,3);
    end
    
    % Elementos que contienen cada nodo
    geom.element2node = element2node(geom.elements);
    % Tabla de conectividad de nodos, bordes, e.t.c
    geom.nodecon2node = node2node(geom.elements);

    % Encuentre los vertices de la malla si se va a usar la adaptacion de malla
    % pasiva de Zinchenco et al. 1997
    geom.vertices = extractvertices(geom);

    % calcule el volumen inicial de la gota
    normalandgeoopt.normal = 1;
    normalandgeoopt.areas = 1;
    normalandgeoopt.vol = 1;
    geomprop = normalandgeo(geom,normalandgeoopt);
    % geom.normalele = geomprop.normalele;
    geom.normal = geomprop.normal;
    geom.dsi = geomprop.dsi;
    geom.ds = geomprop.ds;
    geom.s = geomprop.s;
    geom.vol = geomprop.vol;
    geom.jacmat = geomprop.jacmat;
    geom.volini = geom.vol;
    geom.areaini = geom.s;
    
    % Introduccion de ruido para minimizar efecto de simetria en creacion
    % de la malla.
    
    % Primero se calcula la longitud minima entre los nodos para que el
    % ruido sea una fraccion de esta longitud.
    lmin = zeros(numnodes,1);
    for k = 1:numnodes
       nodesadj = geom.nodecon2node{k};
       lmin(k) = min(normesp(repmat(geom.nodes(k,:),[size(nodesadj,1) 1])...
          - geom.nodes(nodesadj,:)));  
    end

    lmint = min(lmin);

    for k = 1:noiserep
        noisevel = ones(size(geom.nodes))...
            .*(rand(size(geom.nodes))-0.5)*lmint*noiseint;
        noisenormal = repmat(sum(noisevel.*geom.normal,2),[1 3]).*geom.normal;
        noisetan = noisevel - noisenormal;

        geom.nodes = geom.nodes + noisetan;
    end
     
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
    
    normalandgeoopt.normal = 0;
    normalandgeoopt.areas = 1;
    normalandgeoopt.vol = 1;

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
    
    parms.curvopt = parmstemp.curvopt;
    
    normalandgeoopt.normal = 0;
    normalandgeoopt.areas = 1;
    normalandgeoopt.vol = 1;
    
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
    
    normalandgeoopt.normal = 0;
    normalandgeoopt.areas = 1;
    normalandgeoopt.vol = 1;
    
end

% Calculo funciones de forma y demas parametros para el metodo de los 
% elementos finitos que solo depende del estado inicial.

[geom.shapeA, geom.shapeB, geom.refrot] = shapefun(geom);


%% Ciclo principal
abmcount = 0;
      
veladapt = zeros(numnodes,3);

for p = paso:numtimesteps
tic
% if p == 5
%     profile on
% end
% if p == 10
%     profile off
% end
% calcule la distancia minima de adaptacion y el paso de tiempo
    counter = counter + 1;
    disp(['iteracion = ', num2str(p)])
    
%     if p == 4
%        profile on
%     end

    % lmin entre nodos de una misma gota

    lmin = zeros(numnodes,1);
    for k = 1:numnodes
       % extraiga los nodos vecinos a un nodo en la misma gota 
       nodesadj = geom.nodecon2node{k};
       lmin(k) = min(normesp(repmat(geom.nodes(k,:),[size(nodesadj,1) 1])...
          - geom.nodes(nodesadj,:)));  
    end

    % Longitud minima para verificacion de paso de tiempo.
    lmint = min(lmin);
    deltat = lmint^1.5/redfactor;
    parms.lmin = lmin;
    
    %% Para los primeros pasos usamos RK4 para inicializar los puntos de
    %% Adams-Bashforthabmcount = 0;
      if (p <=3 || (opcionsim == 1 && p-paso < 3))
          % primer paso de runge kutta f1
          % invoque el problema de flujo de stokes
          abmcount = abmcount+1;
          [velnode1,geom,parms] = stokesvesicle(geom,parms);
          % invoque la adaptacion de la malla
          % veladapt0 = meshadapt(geom,adaptparms,velnode0);

          % passive (zinchenco et al. 1997)
          velnormal = repmat(sum(velnode1.*geom.normal,2),[1 3]).*geom.normal;
          veladapt = meshadaptgrad(geom,velnormal,veladapt);
          f1 = (velnormal + veladapt);

          abm{abmcount} = f1;
          nodes0 = geom.nodes;

       %% segundo paso de runge kutta f2
          % invoque el problema de flujo de stokes
          geom.nodes = geom.nodes + (1/2)*deltat*f1;
          
          paropt.tipo = 'extended';
          [dummyc,geom.normal,dummyk] = curvparaboloid(geom,paropt);
          
          [velnode2,geom,parms] = stokesvesicle(geom,parms);
          % invoque la adaptacion de la malla
          % veladapt0 = meshadapt(geom,adaptparms,velnode0);

          % passive (zinchenco et al. 1997)
          velnormal = repmat(sum(velnode2.*geom.normal,2),[1 3]).*geom.normal;
          veladapt = meshadaptgrad(geom,velnormal,veladapt);
          f2 = (velnormal + veladapt);

          %% tercer paso de runge kutta f3
          % invoque el problema de flujo de stokes
          geom.nodes = nodes0 + (1/2)*deltat*f2;
          
          paropt.tipo = 'extended';
          [dummyc,geom.normal,dummyk] = curvparaboloid(geom,paropt);
          
          [velnode3,geom,parms] = stokesvesicle(geom,parms);
          % invoque la adaptacion de la malla
          % veladapt0 = meshadapt(geom,adaptparms,velnode0);

          % passive (zinchenco et al. 1997)
          velnormal = repmat(sum(velnode3.*geom.normal,2),[1 3]).*geom.normal;
          veladapt = meshadaptgrad(geom,velnormal,veladapt);
          f3 = (velnormal + veladapt);

       %% cuarto paso de runge kutta f4
          % invoque el problema de flujo de stokes
          geom.nodes = nodes0 + deltat*f3;
          
          paropt.tipo = 'extended';
          [dummyc,geom.normal,dummyk] = curvparaboloid(geom,paropt);
          
          [velnode,geom,parms] = stokesvesicle(geom,parms);
          % invoque la adaptacion de la malla
          % veladapt0 = meshadapt(geom,adaptparms,velnode0);

          % passive (zinchenco et al. 1997)
          velnormal = repmat(sum(velnode.*geom.normal,2),[1 3]).*geom.normal;
          veladapt = meshadaptgrad(geom,velnormal,veladapt);
          f4 = (velnormal + veladapt);

          geom.nodes = nodes0 + deltat*(f1+2*f2+2*f3+f4)/6;
          
          paropt.tipo = 'extended';
          [dummyc,geom.normal,dummyk] = curvparaboloid(geom,paropt);

      else

          % Ahora podemos usar el metodo Predictor-Corrector de
          % Adams-Bashforth-Moulton

          % Paso predictor

          % Calculo de la velocidad en el punto actual

          [velnode1,geom,parms] = stokesvesicle(geom,parms);
          % invoque la adaptacion de la malla
          % veladapt0 = meshadapt(geom,adaptparms,velnode0);

          % passive (zinchenco et al. 1997)
          velnormal = repmat(sum(velnode1.*geom.normal,2),[1 3]).*geom.normal;
          veladapt = meshadaptgrad(geom,velnormal,veladapt);
          f1 = (velnormal + veladapt);

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

          % passive (zinchenco et al. 1997)
          velnormal = repmat(sum(velnode.*geom.normal,2),[1 3]).*geom.normal;
          veladapt = meshadaptgrad(geom,velnormal,veladapt);
          f2 = (velnormal + veladapt);

          abm{5} = f2;

          geom.nodes = nodes0 + deltat*(abm{2}-5*abm{3}+19*abm{4}+9*abm{5})/24;
          
          paropt.tipo = 'extended';
          [dummyc,geom.normal,dummyk] = curvparaboloid(geom,paropt);

          % Actualizaci?n nodos

          abm{1} = abm{2};
          abm{2} = abm{3};
          abm{3} = abm{4};
      end
%   parms.bending.sigma
%% escalaje
    geomprop = normalandgeo(geom,normalandgeoopt);
    geom.dsi = geomprop.dsi;
    geom.ds = geomprop.ds;
    geom.s = geomprop.s;
    geom.vol = geomprop.vol;
    
    errorvol = abs((geom.vol - geom.volini)./geom.volini);

    if errorvol > errorvoltol
           % invoque escalaje
           disp('Realizando Escalaje')
           geom = scaling(geom,optesc,errorvol);
    end

% error de volumen
    errorvol = abs(geom.volini - geom.vol)/geom.volini;
% velocidad normal maxima y tiempo de simulacion
    velcont = max(abs(sum(velnode.*geom.normal,2)));
    geom.tiempo = geom.tiempo + deltat;
    geom.deltat = deltat;

% Visualizacion
    figure(1);
    grafscfld(geom,geom.curv);
    axis equal; view(90,0); xlabel('x1'); ylabel('x2'); zlabel('x3'); colorbar;
    getframe; title('curv');
    figure(2);
    grafscfld(geom,geom.curv);
    axis equal; view(0,0); xlabel('x1'); ylabel('x2'); zlabel('x3'); colorbar;
    getframe; title('curv');
    
  figure(3); grafscfld(geom,geom.lapcurv);
  axis equal; view(90,0); xlabel('x1'); ylabel('x2'); zlabel('x3'); colorbar;
  getframe; title('laplace curv');

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
