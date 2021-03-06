% CALCULO DEL FLUJO DE STOKES PARA UNA VESICULA
% IMPLEMENTADO FLUJO INFINITO, SEMIINFINITO
% IMPLEMENTADO SINGLE Y DOUBLE LAYER
clear;clc; %close all;
%% opciones de carga de archivos
% nombre de archivo a cargar y carpeta
nombreorigen = 'sph ref 5'; %rbc, ellipsoide95
carpetaorigen = '';
iteracion = [];

% nombre de archivo a guardar y carpeta
nombredestino = 'it';
carpetadestino = 'sed_vesicula_hiperelas_grav_inf_maran';
% simulacion nueva desde cero optsim = 0
% continue la simulacion optsim = 1
% simulacion nueva desde archivo de resultados optsim = 2
opcionsim = 0;

% Parametros introduccion ruido para minimizar efecto de la simetria de la
% malla

noiseint = 0.025;
noiserep = 10;

% Algoritmo de flujo de stokes.
ca = 0;
lamda = 1;

% tipo de flujo flow: 'inf'  flow:'semiinf'
flow = 'inf';

% opcion de calculo de la curvatura 1: paraboloid fitting; 2: extended par;
% 3: basado en laplace beltrami
curvopt = 3;

% Coeficientes del modelo de Evans y Skalak
% Coeficiente de resistencia al cambio de area:
% Ka*R_0^2/kappa.
kext = 1e3;
mu = 1;

% gravedad
kb = 1;
g0 = 68.74;

% interaccion electrostatica
ke = 0;
lie = 64.86;
gammaie = 3183.1;

% Coordenadas del centroide de la particula
xc =[0 0 0];

% frecuencia de guardar resultados
outputfreq = 10;

% pasos de tiempo de la simulacion
numtimesteps = 80000;
% Reduccion del p8aso de tiempo calculado automaticamente
redfactor = 5000;

% Sin adaptacion de malla. OJO!
% parametros de adaptacion
% velopt: 1 hidrodinamica velopt:2 normal velopt:3 passive (zinchenko et al.)
% velopt = 3;
% meshadapt lowenberg
% adaptparms.psi = 1;
% adaptparms.lamda = lamda;

% escalaje
errorvoltol = 1e-6;
optesc.maxit = 100;
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

% Coeficiente termino de marangoni
parms.rkmaran = 1/g0;

% Coeficiente adimensional termino bending
parms.rkbend = 1/g0;
parms.kext = kext;
parms.mu = mu;

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

greenfunction = @greeninf;
parms.greenfunction = greenfunction;

dlmod = 1;
parms.dlmod = 'deflaction';

% numero de puntos a usar para integracion polar 4-6-8-12-20
npolar = 4;

[zz,ww] = gausslegabsweights(npolar);
% coordenadas y pesos de los puntos de integracion gauss-Leg 2D
[rmaxh,wwrho,r,xin,etn,ztn] = gausslegintpt(zz,ww);
% parametros de intergacion polar 2D
parms.polarparms.rmaxh = rmaxh;
parms.polarparms.wwrho = wwrho;
parms.polarparms.ww = ww;
parms.polarparms.xin = xin;
parms.polarparms.etn = etn;
parms.polarparms.ztn = ztn;
parms.polarparms.r = r;

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
%     geom.vertices = extractvertices(geom);

    % calcule el volumen inicial de la gota
    normalandgeoopt.normal = 1;
    normalandgeoopt.areas = 1;
    normalandgeoopt.vol = 1;
    geomprop = normalandgeo(geom,normalandgeoopt,1);
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
    
    direcciondestino = ...
        [cd  sbar '..' sbar 'data' sbar carpetadestino sbar nombredestino];
    mkdir([cd  sbar '..' sbar 'data' sbar carpetadestino]);
        
    paso = 1;
    counter = 0;
    geom.tiempo = 0;
    itsaved = 0;
    
    % Geometria de referencia
    geom.ref = geom.nodes;
    geom.dsref = geom.ds;
    
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
    parms = parmstemp;
    counter = 0;
    geom.tiempo = 0;
    itsaved = 0;
    numnodes = size(geom.nodes,1);
    numelements = size(geom.elements,1);
    
    if parms.lamda ~= 1
       geom.W = zeros(numnodes,3); 
       geom.velnodeant = zeros(numnodes,3);
    end
    
    % Geometria de referencia
    geom.ref = geom.nodes;
    geom.dsref = geom.ds;
    
    % volumen reducido inicial
    volredini = geom.volredini;
end

% Calculo funciones de forma y demas parametros para el metodo de los 
% elementos finitos que solo depende del estado inicial.

[geom.shapeA, geom.shapeB, geom.refrot] = shapefun(geom);


%% Ciclo principal
xcant = centroide(geom);
geom.velcentroid = [0 0 0];
geom.xc = xcant;

for p = paso:numtimesteps
% tic
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
    [velnode1,geom,parms] = stokes_wall(geom,parms);
          
    nodes0 = geom.nodes;

    %% segundo paso de runge kutta f2
    % invoque el problema de flujo de stokes
    geom.nodes = nodes0 + (1/2)*deltat*velnode1;

    [velnode,geom,parms] = stokes_wall(geom,parms);

    geom.nodes = nodes0 + deltat*velnode;          

%   parms.bending.sigma
%% escalaje
    normalandgeoopt.normal = 0;
    normalandgeoopt.areas = 1;
    normalandgeoopt.vol = 1;
    geomprop = normalandgeo(geom,normalandgeoopt,1);
    geom.dsi = geomprop.dsi;
    geom.ds = geomprop.ds;
    geom.s = geomprop.s;
    geom.vol = geomprop.vol;
    
    errorvol = abs((geom.vol - geom.volini)./geom.volini);
%     disp(['Error volumen pre-escalaje: ',num2str(errorvol)]);

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

% Visualizacion
    figure(1);
    grafscfld(geom,normesp(geom.rdeltafbend+geom.rdeltafcurv));
    axis equal;
    view(90,0); xlabel('x1'); ylabel('x2'); zlabel('x3'); colorbar;
    title('Normales'); getframe; hold off;
    
    figure(2);
    grafscfld(geom,normesp(geom.rdeltafmaran));
    axis equal;
    view(90,0); xlabel('x1'); ylabel('x2'); zlabel('x3'); colorbar;
    title('Marangoni'); getframe; hold off;
    
    if p > 1
        figure(3); plot(geom.tiempo,geom.velcentroid(3),'*r');hold on;
        title('Sedimentation Rate'); getframe;
    end
%     plot(geom.tiempo,geom.fuerzagrav,'*b'); getframe;
    if ke == 1
        disp(['Fuerzaelest ', num2str(geom.fuerzaelest)]);
    end
    if kb == 1
        disp(['Fuerzagrav: ', num2str(geom.fuerzagrav)]);
    end
    
%     figure(3); plot(geom.tiempo,geom.xc(3),'*k');hold on;
%     title('Posicion Centroide'); getframe;
%         
%     figure(2);
%     grafscfld(geom,normesp(geom.rdeltafmaran));
%     axis equal; view(90,0); xlabel('x1'); ylabel('x2'); zlabel('x3'); colorbar;
%     hold on
%     quiver3(geom.nodes(:,1),geom.nodes(:,2),geom.nodes(:,3),...
%         geom.rdeltafmaran(:,1),geom.rdeltafmaran(:,2),geom.rdeltafmaran(:,3));
%     getframe; title('Marangoni');
%     hold off
    

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
    
    if velcont*deltat < 1e-10
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
