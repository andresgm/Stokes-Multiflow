% CALCULO DEL FLUJO DE STOKES PARA UNA GOTA
% IMPLEMENTADO FLUJO INFINITO, SEMIINFINITO
% IMPLEMENTADO SINGLE Y DOUBLE LAYER
clear;clc;close all;
%% opciones de carga de archivos
    % nombre de archivo a cargar y carpeta
nombreorigen = 'sph ref 3';
carpetaorigen = '';
iteracion = [];
    % nombre de archivo a guardar y carpeta
nombredestino = 'it';
carpetadestino = 'sedimentacion_gota_g0_1_lambda_1_validacion_stokes_inf';
    % simulacion nueva desde cero optsim = 0
    % continue la simulacion optsim = 1
    % simulacion nueva desde archivo de resultados optsim = 2
opcionsim = 0;

% Algoritmo de flujo de stokes.
% capilar o bond
ca = 0;
% lamda
lamda = 1;
% g0: solo aplica para adim = 1. g0 = 1 por defecto
g0 = 1;
% campo electrico
e0 = 0;
% tipo de flujo flow: 'inf'  flow:'semiinf'
flow = 'inf';
% aplica sol cuando hay double layer: 1: 'deflaction' 2:'subsust'
dlmod = 1;
% opcion de calculo de la curvatura 1: paraboloid fitting; 2: best par (extended);
% 3: basado en laplace beltrami
curvopt = 2;
% Adimensionalizacion
adim = 1;
% frecuencia de guardar resultados
outputfreq = 10;

% Banderas de fuerza dif 0: si. 1: no
    % Tension superficial
ka = 1;
    % gravedad
kb = 1;
    % campo electrico
kd = 0;

% numero de gotas
geom.numdrops = 1;
% Coordenadas de los centroides de las gotas
xc =[0 0 20];
% Introduzca el/los radios de la/s gotas
xr=[1];

% pasos de tiempo de la simulacion
numtimesteps = 80000;

redfactor = 70;

% parametros de adaptacion
% velopt: 1 hidrodinamica velopt:2 normal
velopt = 2;
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
parms = conststokesdrop(adim,flow,ca,lamda,g0,e0,ka,kb,kd);
parms.curvopt = curvopt;

if dlmod == 1
    parms.dlmod = 'deflaction';
elseif dlmod == 2
    parms.dlmod = 'subsust';
end

% guarde temporalmente los parametros
parmstemp = parms;
%% procesamiento de la malla
sbar = filesep;
    
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

        % Tabla de elementos singulares a cada nodo
    geom.element2node = element2node(geom.elements);
        % Tabla de conectividad de nodos, bordes, e.t.c  
    geom.nodecon2node = node2node(geom.elements);

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
    geom.areaini = geom.s;
    geom.xcini = centroide(geom);
    
    if parms.lamda ~= 0
        geom.W = zeros(numnodes,3);
        geom.velnodeant = zeros(numnodes,3);
    end
    
         
    if isempty(carpetadestino) == 1
        direccion = [cd  sbar nombredestino num2str(iteracion) '.mat'];   
    else
        direccion = [cd  sbar carpetadestino sbar nombredestino num2str(iteracion) '.mat'];        
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
        direccion = [cd  sbar carpetaorigen sbar nombreorigen num2str(iteracion) '.mat'];
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
    parms.dlmod = parmstemp.dlmod;
    
elseif opcionsim == 2
    % cargue desde resultados y realice una nueva simulacion
    % cargue desde resultados y continue la simulacion
    if isempty(carpetaorigen) == 1
        direccion = [cd  sbar nombreorigen num2str(iteracion) '.mat'];        
        
    else
        direccion = [cd  sbar carpetaorigen sbar nombreorigen num2str(iteracion) '.mat'];        
    end
    
    load(direccion);
    
    direcciondestino = [cd  sbar carpetadestino sbar nombredestino];
    mkdir([cd sbar,carpetadestino]);
        
    paso = 1;    
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


%% evolucion mediante runge kutta de 4to orden
xcant = centroide(geom);
geom.velcentroid = [0 0 0];
geom.xc = xcant;
for p = paso:numtimesteps
tic
% calcule la distancia minima de adaptacion y el paso de tiempo
    counter = counter + 1;
    disp(['iteracion = ', num2str(p)])
   
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
    parms.lmin = lmin;
    
% Para los primeros pasos usamos RK4 para inicializar los puntos de
% Adams-Bashforth-Moulton
   if p <=3
       % primer paso de runge kutta f1
       % invoque el problema de flujo de stokes
       [velnode1,geom,parms] = stokesdrop(geom,parms);
       % invoque la adaptacion de la malla
       % veladapt0 = meshadapt(geom,adaptparms,velnode0);

       if velopt == 1
           % hidrodinamica + adaptacion
           veladapt1 = meshadapt2(geom,parms);
           temporal = repmat(geom.velcentroid,[numnodes 1]);
    %        veltan = temporal - repmat(sum(temporal.*geom.normal,2),[1 3]).*geom.normal;
           veltan = 0;
           f1 = deltat.*(velnode1 + veladapt1 + veltan);
       elseif velopt == 2
           % normal + adaptacion
           temporal = repmat(geom.velcentroid,[numnodes 1]);
           veltan = temporal - repmat(sum(temporal.*geom.normal,2),[1 3]).*geom.normal;
    %        veltan = 0;
    %        veladapt = 0;
           velnormal1 = repmat(sum(velnode1.*geom.normal,2),[1 3]).*geom.normal;
           veladapt1 = meshadapt2(geom,parms);
           f1 = deltat.*(velnormal1 + veladapt1 + veltan);
       end

       abm{counter} = f1;
       nodes0 = geom.nodes;
       geom.nodes = geom.nodes + (1/2).*f1;

    %% segundo paso de runge kutta f2
       % invoque el problema de flujo de stokes
       [velnode2,geom,parms] = stokesdrop(geom,parms);
       % invoque la adaptacion de la malla
    %    veladapt0 = meshadapt(geom,adaptparms,velnode0);

       if velopt == 1
           % hidrodinamica + adaptacion
           veladapt2 = meshadapt2(geom,parms);
           temporal = repmat(geom.velcentroid,[numnodes 1]);
    %        veltan = temporal - repmat(sum(temporal.*geom.normal,2),[1 3]).*geom.normal;
           veltan = 0;
           f2 = deltat.*(velnode2 + veladapt2 + veltan);
       elseif velopt == 2
           % normal + adaptacion
           temporal = repmat(geom.velcentroid,[numnodes 1]);
           veltan = temporal - repmat(sum(temporal.*geom.normal,2),[1 3]).*geom.normal;
    %        veltan = 0;
    %        veladapt = 0;
           velnormal2 = repmat(sum(velnode2.*geom.normal,2),[1 3]).*geom.normal;
           veladapt2 = meshadapt2(geom,parms);
           f2 = deltat.*(velnormal2 + veladapt2 + veltan);
       end

       geom.nodes = nodes0 + (1/2).*f2;

    %% tercer paso de runge kutta f3
       % invoque el problema de flujo de stokes
       [velnode3,geom,parms] = stokesdrop(geom,parms);
       % invoque la adaptacion de la malla
    %    veladapt0 = meshadapt(geom,adaptparms,velnode0);

       if velopt == 1
           % hidrodinamica + adaptacion
           veladapt3 = meshadapt2(geom,parms);
           temporal = repmat(geom.velcentroid,[numnodes 1]);
    %        veltan = temporal - repmat(sum(temporal.*geom.normal,2),[1 3]).*geom.normal;
           veltan = 0;
           f3 = deltat.*(velnode3 + veladapt3 + veltan);
       elseif velopt == 2
           % normal + adaptacion
           temporal = repmat(geom.velcentroid,[numnodes 1]);
           veltan = temporal - repmat(sum(temporal.*geom.normal,2),[1 3]).*geom.normal;
    %        veltan = 0;
    %        veladapt = 0;
           velnormal3 = repmat(sum(velnode3.*geom.normal,2),[1 3]).*geom.normal;
           veladapt3 = meshadapt2(geom,parms);
           f3 = deltat.*(velnormal3 + veladapt3 + veltan);
       end

       geom.nodes = nodes0 + f3;

    %% cuarto paso de runge kutta f4
       % invoque el problema de flujo de stokes
       [velnode,geom,parms] = stokesdrop(geom,parms);
       % invoque la adaptacion de la malla
    %    veladapt0 = meshadapt(geom,adaptparms,velnode0);

       if velopt == 1
           % hidrodinamica + adaptacion
           veladapt = meshadapt2(geom,parms);
           temporal = repmat(geom.velcentroid,[numnodes 1]);
    %        veltan = temporal - repmat(sum(temporal.*geom.normal,2),[1 3]).*geom.normal;
           veltan = 0;
           f4 = deltat.*(velnode + veladapt + veltan);
       elseif velopt == 2
           % normal + adaptacion
           temporal = repmat(geom.velcentroid,[numnodes 1]);
           veltan = temporal - repmat(sum(temporal.*geom.normal,2),[1 3]).*geom.normal;
    %        veltan = 0;
    %        veladapt = 0;
           velnormal = repmat(sum(velnode.*geom.normal,2),[1 3]).*geom.normal;
           veladapt = meshadapt2(geom,parms);
           f4 = deltat.*(velnormal + veladapt + veltan);
       end

       geom.nodes = nodes0 + (f1+2*f2+2*f3+f4)/6;

   else
       
       % Ahora podemos usar el m?todo Predictor-Corrector de
       % Adams-Bashforth-Moulton
       
       % Paso predictor
       
       % Calculo de la velocidad en el punto actual
       
       [velnode1,geom,parms] = stokesdrop(geom,parms);
       % invoque la adaptacion de la malla
       % veladapt0 = meshadapt(geom,adaptparms,velnode0);

       if velopt == 1
           % hidrodinamica + adaptacion
           veladapt1 = meshadapt2(geom,parms);
           temporal = repmat(geom.velcentroid,[numnodes 1]);
    %        veltan = temporal - repmat(sum(temporal.*geom.normal,2),[1 3]).*geom.normal;
           veltan = 0;
           f1 = deltat.*(velnode1 + veladapt1 + veltan);
       elseif velopt == 2
           % normal + adaptacion
           temporal = repmat(geom.velcentroid,[numnodes 1]);
           veltan = temporal - repmat(sum(temporal.*geom.normal,2),[1 3]).*geom.normal;
    %        veltan = 0;
    %        veladapt = 0;
           velnormal1 = repmat(sum(velnode1.*geom.normal,2),[1 3]).*geom.normal;
           veladapt1 = meshadapt2(geom,parms);
           f1 = deltat.*(velnormal1 + veladapt1 + veltan);
       end

       abm{4} = f1;
       nodes0 = geom.nodes;
       geom.nodes = geom.nodes + (-9*abm{1}+37*abm{2}-59*abm{3}+55*abm{4})/24;
       
       % Paso corrector
       
       % Calculo de la velocidad en el punto siguiente con la prediccion
       
       [velnode,geom,parms] = stokesdrop(geom,parms);
       % invoque la adaptacion de la malla
       % veladapt0 = meshadapt(geom,adaptparms,velnode0);

       if velopt == 1
           % hidrodinamica + adaptacion
           veladapt = meshadapt2(geom,parms);
           temporal = repmat(geom.velcentroid,[numnodes 1]);
    %        veltan = temporal - repmat(sum(temporal.*geom.normal,2),[1 3]).*geom.normal;
           veltan = 0;
           f2 = deltat.*(velnode + veladapt + veltan);
       elseif velopt == 2
           % normal + adaptacion
           temporal = repmat(geom.velcentroid,[numnodes 1]);
           veltan = temporal - repmat(sum(temporal.*geom.normal,2),[1 3]).*geom.normal;
    %        veltan = 0;
    %        veladapt = 0;
           velnormal = repmat(sum(velnode.*geom.normal,2),[1 3]).*geom.normal;
           veladapt = meshadapt2(geom,parms);
           f2 = deltat.*(velnormal + veladapt + veltan);
       end

       abm{5} = f2;
       
       geom.nodes = nodes0 + (abm{2}-5*abm{3}+19*abm{4}+9*abm{5})/24;
       
       % Actualizaci?n nodos
       
       abm{1} = abm{2};
       abm{2} = abm{3};
       abm{3} = abm{4};
       
   
   end

%% escalaje
    geomprop = normalandgeo(geom,normalandgeoopt);
    geom.normalele = geomprop.normalele;
    geom.normal = geomprop.normal;
    geom.dsi = geomprop.dsi;
    geom.ds = geomprop.ds;
    geom.s = geomprop.s;
    geom.vol = geomprop.vol;
    geom.jacmat = geomprop.jacmat;
    
    errorvol = abs((geom.vol - geom.volini)./geom.volini);

    for r=1:geom.numdrops
       if errorvol(r) > errorvoltol
           % invoque escalaje
           disp('Realizando Escalaje')
           geom = escaling(geom,optesc,r);
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
    
% grafique la geometria
    figure(1);
    grafscfld(geom,geom.curv);
    axis equal; view(90,0); xlabel('x1'); ylabel('x2'); zlabel('x3'); colorbar;
    getframe; title('curv');
% % grafique la velocidad del centroide
%     figure(2); plot(geom.tiempo,normesp(geom.velcentroid),'*');hold on; title('velcentroid');
% % grafique la velocidad normal
    figure(3); plot(geom.tiempo,velcont,'*');hold on; title('velnormal');
% % grafique error de volumen
%     figure(4); plot(geom.tiempo,errorvol,'*');hold on;title('errorvol');
% grafique la velocidad del centroide en x3
    figure(4); plot(geom.tiempo,geom.velcentroid(:,3),'*');hold on; title('vc in x3');
% grafique la posicion del centroide en x3
    figure(5); plot(geom.tiempo,geom.xc(:,3),'*');hold on; title('xc x3');   
    
% guarde resultados
    if counter == outputfreq
        itsaved = itsaved + 1;
        counter = 0;
        nombrearchivo = [direcciondestino num2str(itsaved), '.mat'];
        save(nombrearchivo,'geom','velnode','parms','adim');
    end
    disp(carpetadestino)
    disp(['deltat: ', num2str(geom.deltat)])
    
toc
end
