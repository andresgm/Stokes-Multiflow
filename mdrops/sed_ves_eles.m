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
carpetadestino = 'sed_vesicula_g01_fig4';
% simulacion nueva desde cero optsim = 0
% continue la simulacion optsim = 1
% simulacion nueva desde archivo de resultados optsim = 2
opcionsim = 0;

% Algoritmo de flujo de stokes.
ca = 0;
lamda = 1;

% tipo de flujo flow: 'inf'  flow:'semiinf'
flow = 'semiinf';

% opcion de calculo de la curvatura 1: paraboloid fitting; 2: extended par;
% 3: basado en laplace beltrami
curvopt = 3;

% Constantes del modelo de Evans y Rawicz
% Constante c del modelo
c = 0.1;
% Coeficiente de rigides al doblamiento en KbT (kappa)
kbar = 25;
% Coeficiente adimensional de resistencia al cambio de area:
% Ka*R_0^2/kappa.
kext = 7.58e8;

% Banderas de fuerza dif 0: si. 1: no
% gravedad
kb = 1;
g0 = 1;

% interaccion electrostatica
ke = 1;
lie = 33.33;
gammaie = 150.71;

% Coordenadas de los centroides de las gotas
xc =[0 0 5];

% frecuencia de guardar resultados
outputfreq = 10;

% pasos de tiempo de la simulacion
numtimesteps = 80000;
% Reduccion del paso de tiempo calculado automaticamente
redfactor = 100;

% parametros de adaptacion
% velopt: 1 hidrodinamica velopt:2 normal velopt:3 passive (zinchenko et al.)
velopt = 3;

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
parms.rkextf = 2*ca/g0;
parms.rksl = 2;
parms.rkdl = 2*(lamda - 1)/(lamda + 1);
parms.lamda = lamda;
parms.ca = ca;
parms.g0 = g0;   

% Coeficiente termino de curvatura
parms.rkcurv = 1/g0;
% Coeficiente adimensional termino bending
% bending
parms.rkbend = 1/g0;
parms.bending.c = c;
parms.bending.kbar = kbar;
parms.bending.kext = kext;
parms.bending.sigma = 0;

if kb == 1
    % gravedad
    parms.rkgrav = 1;
else
    parms.rkgrav = 0;
end

if ke == 1
   % Interaccion electrostatica
   parms.rkelestat = gammaie/g0;
   parms.elestat.l = lie;
else
    parms.rkelestat = 0;
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
    % Numero de Elementos y numero de Nodos
    geom.nodes = Nodes;
    geom.elements = Elements;
    geom.numnodes = size(geom.nodes,1);
    geom.numelements = size(Elements,1);
    numnodes = geom.numnodes;
    numelements = geom.numelements;
    geom = set_centroid(geom,xc);
    
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
    
    direcciondestino = ...
        [cd  sbar '..' sbar 'data' sbar carpetadestino sbar nombredestino];
    mkdir([cd  sbar '..' sbar 'data' sbar carpetadestino]);
        
    paso = 1;
    counter = 0;
    geom.tiempo = 0;
    itsaved = 0;
    
    % volumen reducido inicial
    volredini = 6*sqrt(pi)*geom.vol/geom.s^(3/2);
    geom.volredini = volredini;
    disp(['Volumen reducido incial: ',num2str(volredini)]);

elseif opcionsim == 1
    % cargue desde resultados y continue la simulacion
    carpetadestino = carpetaorigen;
    nombredestino = nombreorigen;
    if isempty(carpetaorigen) == 1
        error('Must specify folder to read previous result.');        
    else
        direccion = ...
         [cd  sbar '..' sbar 'data' sbar carpetaorigen sbar nombreorigen...
         num2str(iteracion) '.mat'];
    end
    
    load(direccion);
    
    paso = iteracion + 1;
    counter = 0;
    itsaved = iteracion;

    direcciondestino = ...
        [cd  sbar '..' sbar 'data' sbar carpetadestino sbar nombredestino];
    mkdir([cd  sbar '..' sbar 'data' sbar carpetadestino]);
    numnodes = size(geom.nodes,1);
    numelements = size(geom.elements,1);    
    
    parms.curvopt = parmstemp.curvopt;
    % volumen reducido inicial
    volredini = geom.volredini;
    
elseif opcionsim == 2
    % cargue desde resultados y realice una nueva simulacion
    % cargue desde resultados y continue la simulacion
    if isempty(carpetaorigen) == 1
        error('Must specify folder to read previous result.');        
    else
        direccion = ...
         [cd  sbar '..' sbar 'data' sbar carpetaorigen sbar nombreorigen...
         num2str(iteracion) '.mat'];
    end
    
    load(direccion);
    
    direcciondestino = ...
        [cd  sbar '..' sbar 'data' sbar carpetadestino sbar nombredestino];
    mkdir([cd  sbar '..' sbar 'data' sbar carpetadestino]);
        
    paso = 1;    
    parmstemp.bending.sigma = parms.bending.sigma;
    parms = parmstemp;
    counter = 0;
    geom.tiempo = 0;
    itsaved = 0;
    numnodes = size(geom.nodes,1);
    numelements = size(geom.elements,1);
    
    % volumen reducido inicial
    volredini = geom.volredini;
end


%% Ciclo principal
xcant = centroide(geom);
geom.velcentroid = [0 0 0];
geom.xc = xcant;
      
veladapt = zeros(numnodes,3);

for p = paso:numtimesteps
%tic
% if p == 5
%     profile on
% end
% if p == 10
%     profile off
% end

    counter = counter + 1;
    disp('************************')
    disp(['iteracion = ', num2str(p)])

% calcule la distancia minima de adaptacion y el paso de tiempo

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
    
    % Usamos esquema rk2 para la integracion numerica.
    
    %% primer paso de runge kutta f1
    [velnode1,geom,parms] = stokesvesicle_sed_eles(geom,parms);
    
    % passive (zinchenco et al. 1997)
    velnormal = repmat(sum(velnode1.*geom.normal,2),[1 3]).*geom.normal;
    veladapt = meshadaptgrad(geom,velnormal,veladapt);
    f1 = (velnormal + veladapt);
          
    nodes0 = geom.nodes;
    
    geom.nodes = nodes0 + (1/2)*deltat*velnode1;
    
    %% segundo paso de runge kutta f2
    % invoque el problema de flujo de stokes

    [velnode,geom,parms] = stokesvesicle_sed_eles(geom,parms);
    
    % passive (zinchenco et al. 1997)
    velnormal = repmat(sum(velnode.*geom.normal,2),[1 3]).*geom.normal;
    veladapt = meshadaptgrad(geom,velnormal,veladapt);
    f2 = (velnormal + veladapt);
       
    geom.nodes = nodes0 + deltat*f2;

%   parms.bending.sigma
%% escalaje
    normalandgeoopt.normal = 0;
    normalandgeoopt.areas = 1;
    normalandgeoopt.vol = 1;
    geomprop = normalandgeo(geom,normalandgeoopt);
    geom.dsi = geomprop.dsi;
    geom.ds = geomprop.ds;
    geom.s = geomprop.s;
    geom.vol = geomprop.vol;
    
    errorvol = abs((geom.vol - geom.volini)./geom.volini);

    if errorvol > errorvoltol
       % invoque escalaje
       geom = scaling(geom,optesc,errorvol);
    end

% error de volumen
%     errorvol = abs(geom.volini - geom.vol)/geom.volini;
%     disp(['Error volumen post-escalaje: ',num2str(errorvol)]);
% velocidad normal maxima y tiempo de simulacion

    velcont = max(abs(sum(velnode.*geom.normal,2)));
    disp(['Velocidad normal maxima: ', num2str(velcont)]);
    geom.tiempo = geom.tiempo + deltat;
    geom.deltat = deltat;
    
% calculo del centroide de la gota y velocidad del centroide
    xcant = geom.xc;
    geom.xc = centroide(geom);
    geom.velcentroid = (geom.xc - xcant)./deltat;
    disp(['Posicion centroide: ', num2str([geom.xc(1),geom.xc(2),geom.xc(3)])]);
    disp(['Velocidad centroide: ', num2str(geom.velcentroid)]);
    % todo: revisar este display
    
% Visualizacion
    figure(1);
    grafscfld(geom,geom.rdeltafnorm);
    axis equal;
    view(90,0); xlabel('x1'); ylabel('x2'); zlabel('x3'); colorbar;
    title('Tension normal'); getframe; hold off;
    
%     figure(2); plot(geom.tiempo,parms.bending.sigma,'*'); hold on;title('sigma');
    
    if ke == 1
        disp(['Fuerzaelest ', num2str(geom.fuerzaelest)]);
    end
    if kb == 1
        disp(['Fuerzagrav: ', num2str(geom.fuerzagrav)]);
    end
    
% guarde resultados
    if counter == outputfreq
        itsaved = itsaved + 1;
        counter = 0;
        nombrearchivo = [direcciondestino num2str(itsaved), '.mat'];
        save(nombrearchivo,'geom','velnode','parms');
    end
%     disp(carpetadestino)
    disp(['deltat: ', num2str(geom.deltat)]);
    disp(['tiempo: ', num2str(geom.tiempo)]);
    volred = 6*sqrt(pi)*geom.vol/geom.s^(3/2);
%     disp(['Volumen reducido: ',num2str(volred)]);
    errorvolred = abs(volred-volredini)/volredini;
    disp(['Error volumen reducido: ',num2str(errorvolred)]);
    
    if velcont*deltat < 1e-6
        disp('Convergencia a estado estacionario');
        itsaved = itsaved + 1;
        nombrearchivo = [direcciondestino num2str(itsaved), '.mat'];
        save(nombrearchivo,'geom','velnode','parms');
        break;
    end
    
%     if p == 14
%        profile viewer
%     end
    
% toc
end
