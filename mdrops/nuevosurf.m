% CALCULO DEL FLUJO DE STOKES PARA UNA GOTA
% IMPLEMENTADO FLUJO INFINITO, SEMIINFINITO
% IMPLEMENTADO SINGLE Y DOUBLE LAYER
clear;clc; close all
%% opciones de carga de archivos
    % nombre de archivo a cargar y carpeta
nombreorigen = 'sph ref 3';
carpetaorigen = '';
iteracion = [];
    % nombre de archivo a guardar y carpeta
nombredestino = 'it';
carpetadestino = 'sedimentacion_gota_g0_1_lambda_1_validacion_stokes_surf_inf';
    % simulacion nueva desde cero optsim = 0
    % continue la simulacion optsim = 1
    % simulacion nueva desde archivo de resultados optsim = 2
opcionsim = 0;

% Algoritmo de flujo de stokes.
% capilar o bond
ca = 0;
% lamda
lamda = 1;
% tipo de flujo flow: 'inf'  flow:'semiinf'
flow = 'inf';
% aplica sol cuando hay double layer: 1: 'deflaction' 2:'subsust'
dlmod = 1;
% opcion de calculo de la curvatura 1: paraboloid fitting; 2: extended par;
% 3: basado en laplace beltrami
curvopt = 1;
% Adimensionalizacion
adim = 2;
% frecuencia de guardar resultados
outputfreq = 10;

% Banderas de fuerza dif 0: si. 1: no
    % Tension superficial
ka = 1;
    % gravedad
kb = 1;
   % Si hay gravedad se debe definir el n?mero de Bond y admin debe ser 2!
Bo = 1;
    % marangoni
kc = 1;
    
% sulfactantes si kc = 1. 
% maranmodel = 1(lineal) definir beta y pe(alpha)
% maranmodel = 2(logaritmico) definir x, e y pe(alpha)
 maranmodel = 2;
 e = 0.20;
 % Concentracion 
 x = 0.975;
 % Parametro \Alpha =  \simga_0 R_0/\mu D_s
 alpha = 100;
 % constante para modelo lineal
 beta = 0.2;

% numero de gotas
geom.numdrops = 1;
% Coordenadas de los centroides de las gotas
xc =[0 0 20];
% Introduzca el/los radios de la/s gotas
xr = 1;

% pasos de tiempo de la simulacion
numtimesteps = 80000;

redfactor = 50;

% parametros de adaptacion
% velopt: 1 hidrodinamica velopt:2 normal velopt:3 passive (zinchenko et al.)
% meshadapt lowenberg
adaptparms.psi = 1;
adaptparms.lamda = lamda;
% escalaje
errorvoltol = 1e-6;
optesc.maxit = 15000;
optesc.kp = 20;
optesc.deltate = 0.01;
optesc.tolerrorvol = errorvoltol;

% numero de puntos a usar para integracion polar 4-6-8-12-20
npolar = 4;

% parametros de surfactantes
% opciones de evolucion de la interfase
% surfopt.opt = 1 los nodos solo se mueven con la velocidad normal
% surfopt.opt = 2 los nodos se mueven con la velocidad normal y adaptacion
% surfopt.opt = 3 los nodos se mueven con la velocidad hidrodinamica
% surfopt.opt = 4 los nodos se mueven con la velocidad hidrodinamica y adaptacion
% surfopt.opt = 5 passive (zinchenko et al.)
surfopt.opt = 5;
% parametros de tiempo para integracion de surfactantes
% theta = 0 Euler Explicito
% theta = 0.5 semi implicito
% theta = 1 full implicito
theta = 0.5;

%% procesamiento de parametros
if adim == 1
    warning('Esta dimensionalizacion permite solo flujo infinito')
    greenfunction = @greeninf;
    parms.greenfunction = greenfunction;
    parms.flow = flow;
    parms.curvopt = curvopt;
       
   % adimensionalizacion de bazhlekov
    parms.rkextf = 2/(lamda + 1);
    parms.rksl = 2/(lamda + 1);
    parms.rkdl = 2*(lamda - 1)/(lamda + 1);
    
    if ca == 0
       error('Para esta adimensionalizacion Ca debe ser un n?mero positivo > 0')
    else
       parms.ca = ca;
    end
    
    parms.lamda = lamda;
    % Peclet para la evolucion de surfactantes
    % pe = ca * ALPHA (Bazhekov)
    s0seq = (1 + e*log(1-x))^(-1);
    % para modelo no lineal    
    if maranmodel == 1
        % lineal
        pe = alpha*ca; 
    elseif maranmodel == 2
        % logaritmico
        pe = ca*alpha*s0seq;
%         pe = alpha
    end
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
    % adimensionalizacion de andres gonzalez basado en la velocidad
    % caracter?stica de sedimentacion U_0=\Delata \rho g R_0^2/\mu(1+\lambda)
    
    parms.flow = flow;
    greenfunction = @greeninf;
    parms.greenfunction = greenfunction;
    parms.curvopt = curvopt;
    % adimensionalizacion del single layer
    parms.rkextf = 2*ca/Bo;
    parms.rksl = 2;
    parms.rkdl = 2*(lamda - 1)/(lamda + 1);

    parms.Bo = Bo;   
    parms.lamda = lamda;
    parms.ca = ca;
    % parametros de simulacion
    if ka == 1 
        % curvatura constante
        parms.rkcurv = 1/Bo;
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
       pe = Bo*alpha/(1-lamda);
        parms.maran.rkmaran = 1/Bo;
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
      
end

if dlmod == 1
    parms.dlmod = 'deflaction';
elseif dlmod == 2
    parms.dlmod = 'subsust';
end

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
        % Tabla de elementos singulares a cada nodo -sparse-
    % geom.e2nsparse = ele2nodesp(geom.elements);
    
    % Encuentre los vertices de la malla si se va a usar la adaptacion de malla
    % pasiva de Zinchenco et al. 1997
    if surfopt.opt == 5
       geom.vertices = extractvertices(geom);
    end

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
    
    % Campo inicial de concentracion
    flds.gamma = ones(numnodes,1);
    geom.gammatotori = inttrapecioa(geom.dsi,flds.gamma)./geom.s;
%     flds.gamma = abs(geom.nodes(:,3));

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
if surfopt.opt == 5
   veladapt = zeros(numnodes,3);
end
for p = 1:numtimesteps
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
       [velnode1,geom] = stokessurf(geom,parms,flds);
       % invoque la adaptacion de la malla
       % veladapt0 = meshadapt(geom,adaptparms,velnode0);

       % factor k1 de runge kutta y calculo de las matrices de concentracion
    if surfopt.opt == 1
       % normal
       velnormal0 = repmat(sum(velnode1.*geom.normal,2),[1 3]).*geom.normal;
       if parms.maran.rkmaran ~= 0
       ajimat = surfactants(geom,velnode1,zeros(numnodes,3),pe,surfopt);
       end
       f1 = velnormal0;
    elseif surfopt.opt == 2
       % normal + adaptacion
       temporal = repmat(geom.velcentroid,[numnodes 1]);
       veltan = temporal - repmat(sum(temporal.*geom.normal,2),[1 3]).*geom.normal;
       veladapt0 = meshadapt2(geom,parms);
       % veladapt0 = veladapt0 - repmat(sum(veladapt0.*geom.normal,2),[1 3]).*geom.normal; 
       velnormal0 = repmat(sum(velnode1.*geom.normal,2),[1 3]).*geom.normal;
       if parms.maran.rkmaran ~= 0
         ajimat = surfactants(geom,velnode1,veladapt0,pe,surfopt);
       end
       f1 = (velnormal0 + veladapt0 + veltan);
    elseif surfopt.opt == 3
       % hidrodinamica
       if parms.maran.rkmaran ~= 0
       ajimat = surfactants(geom,velnode1,zeros(numnodes,3),pe,surfopt);
       end
       f1 = velnode1;
    elseif surfopt.opt == 4
       % hidrodinamica + adaptacion      
       veladapt0 = meshadapt2(geom,parms);
       if parms.maran.rkmaran ~= 0
       ajimat = surfactants(geom,velnode1,veladapt0,pe,surfopt);
       end
       f1 = (velnode1 + veladapt0);
    elseif surfopt.opt == 5
           % passive (zinchenco et al. 1997)
           velnormal = repmat(sum(velnode1.*geom.normal,2),[1 3]).*geom.normal;
           veladapt = meshadaptgrad(geom,velnormal,veladapt);
       if parms.maran.rkmaran ~= 0
         ajimat = surfactants(geom,velnode1,veladapt,pe,surfopt);
       end
       f1 = (velnormal + veladapt);
    end
   
    if parms.maran.rkmaran ~= 0
        % evolucion del sulfactante en t + dt/2
        % propuesto
        %     gammaori = flds.gamma;
        flds.gamma = thetamethod(ajimat,deltat/2,theta,flds.gamma);

        % escale la concentracion
        gammatotnew = inttrapecioa(geom.dsi,flds.gamma)./geom.s;
        gammasc = geom.gammatotori/gammatotnew;
        flds.gamma = flds.gamma.*gammasc;
        geom.gammatot = inttrapecioa(geom.dsi,flds.gamma)./geom.s;
        geom.gammatot
    end

       abm{counter} = f1;
       nodes0 = geom.nodes;
       
    %% segundo paso de runge kutta f2
       % invoque el problema de flujo de stokes
       geom.nodes = geom.nodes + (1/2)*deltat*f1;
       [velnode2,geom] = stokessurf(geom,parms,flds);
       % invoque la adaptacion de la malla
    %    veladapt0 = meshadapt(geom,adaptparms,velnode0);

       % factor k1 de runge kutta y calculo de las matrices de concentracion
    if surfopt.opt == 1
       % normal
       velnormal0 = repmat(sum(velnode2.*geom.normal,2),[1 3]).*geom.normal;
       if parms.maran.rkmaran ~= 0
       ajimat = surfactants(geom,velnode2,zeros(numnodes,3),pe,surfopt);
       end
       f2 = velnormal0;
    elseif surfopt.opt == 2
       % normal + adaptacion
       temporal = repmat(geom.velcentroid,[numnodes 1]);
       veltan = temporal - repmat(sum(temporal.*geom.normal,2),[1 3]).*geom.normal;
       veladapt0 = meshadapt2(geom,parms);
       % veladapt0 = veladapt0 - repmat(sum(veladapt0.*geom.normal,2),[1 3]).*geom.normal; 
       velnormal0 = repmat(sum(velnode2.*geom.normal,2),[1 3]).*geom.normal;
       if parms.maran.rkmaran ~= 0
       ajimat = surfactants(geom,velnode2,veladapt0,pe,surfopt);
       end
       f2 = (velnormal0 + veladapt0 + veltan);
    elseif surfopt.opt == 3
       % hidrodinamica
       if parms.maran.rkmaran ~= 0
       ajimat = surfactants(geom,velnode2,zeros(numnodes,3),pe,surfopt);
       end
       f2 = velnode0;
    elseif surfopt.opt == 4
       % hidrodinamica + adaptacion      
       veladapt0 = meshadapt2(geom,parms);
       if parms.maran.rkmaran ~= 0
       ajimat = surfactants(geom,velnode2,veladapt0,pe,surfopt);
       end
       f2 = (velnode2 + veladapt0);
    elseif surfopt.opt == 5
           % passive (zinchenco et al. 1997)
           velnormal = repmat(sum(velnode2.*geom.normal,2),[1 3]).*geom.normal;
           veladapt = meshadaptgrad(geom,velnormal,veladapt);
       if parms.maran.rkmaran ~= 0
         ajimat = surfactants(geom,velnode2,veladapt,pe,surfopt);
       end
       f2 = (velnormal + veladapt);
    end
   
    if parms.maran.rkmaran ~= 0
        % evolucion del sulfactante en t + dt/2
        % propuesto
        %     gammaori = flds.gamma;
        flds.gamma = thetamethod(ajimat,deltat/2,theta,flds.gamma);

        % escale la concentracion
        gammatotnew = inttrapecioa(geom.dsi,flds.gamma)./geom.s;
        gammasc = geom.gammatotori/gammatotnew;
        flds.gamma = flds.gamma.*gammasc;
        geom.gammatot = inttrapecioa(geom.dsi,flds.gamma)./geom.s;
        geom.gammatot
    end

    %% tercer paso de runge kutta f3
       % invoque el problema de flujo de stokes
       geom.nodes = nodes0 + (1/2)*deltat*f2;
       [velnode3,geom] = stokessurf(geom,parms,flds);
       % invoque la adaptacion de la malla
    %    veladapt0 = meshadapt(geom,adaptparms,velnode0);

       % factor k1 de runge kutta y calculo de las matrices de concentracion
    if surfopt.opt == 1
       % normal
       velnormal0 = repmat(sum(velnode3.*geom.normal,2),[1 3]).*geom.normal;
       if parms.maran.rkmaran ~= 0
       ajimat = surfactants(geom,velnode3,zeros(numnodes,3),pe,surfopt);
       end
       f3 = velnormal0;
    elseif surfopt.opt == 2
       % normal + adaptacion
       temporal = repmat(geom.velcentroid,[numnodes 1]);
       veltan = temporal - repmat(sum(temporal.*geom.normal,2),[1 3]).*geom.normal;
       veladapt0 = meshadapt2(geom,parms);
       % veladapt0 = veladapt0 - repmat(sum(veladapt0.*geom.normal,2),[1 3]).*geom.normal; 
       velnormal0 = repmat(sum(velnode3.*geom.normal,2),[1 3]).*geom.normal;
       if parms.maran.rkmaran ~= 0
       ajimat = surfactants(geom,velnode3,veladapt0,pe,surfopt);
       end
       f3 = (velnormal0 + veladapt0 + veltan);
    elseif surfopt.opt == 3
       % hidrodinamica
       if parms.maran.rkmaran ~= 0
       ajimat = surfactants(geom,velnode3,zeros(numnodes,3),pe,surfopt);
       end
       f3 = velnode3;
    elseif surfopt.opt == 4
       % hidrodinamica + adaptacion      
       veladapt0 = meshadapt2(geom,parms);
       if parms.maran.rkmaran ~= 0
       ajimat = surfactants(geom,velnode3,veladapt0,pe,surfopt);
       end
       f3 = (velnode3 + veladapt0);
    elseif surfopt.opt == 5
           % passive (zinchenco et al. 1997)
           velnormal = repmat(sum(velnode3.*geom.normal,2),[1 3]).*geom.normal;
           veladapt = meshadaptgrad(geom,velnormal,veladapt);
       if parms.maran.rkmaran ~= 0
         ajimat = surfactants(geom,velnode3,veladapt,pe,surfopt);
       end
       f3 = (velnormal + veladapt);
    end
   
    if parms.maran.rkmaran ~= 0
        % evolucion del sulfactante en t + dt/2
        % propuesto
        %     gammaori = flds.gamma;
        flds.gamma = thetamethod(ajimat,deltat/2,theta,flds.gamma);

        % escale la concentracion
        gammatotnew = inttrapecioa(geom.dsi,flds.gamma)./geom.s;
        gammasc = geom.gammatotori/gammatotnew;
        flds.gamma = flds.gamma.*gammasc;
        geom.gammatot = inttrapecioa(geom.dsi,flds.gamma)./geom.s;
        geom.gammatot
    end

    %% cuarto paso de runge kutta f4
       % invoque el problema de flujo de stokes
       geom.nodes = nodes0 + deltat*f3;
       [velnode,geom] = stokessurf(geom,parms,flds);
       % invoque la adaptacion de la malla
    %    veladapt0 = meshadapt(geom,adaptparms,velnode0);

       % factor k1 de runge kutta y calculo de las matrices de concentracion
    if surfopt.opt == 1
       % normal
       velnormal0 = repmat(sum(velnode.*geom.normal,2),[1 3]).*geom.normal;
       if parms.maran.rkmaran ~= 0
       ajimat = surfactants(geom,velnode,zeros(numnodes,3),pe,surfopt);
       end
       f4 = velnormal0;
    elseif surfopt.opt == 2
       % normal + adaptacion
       temporal = repmat(geom.velcentroid,[numnodes 1]);
       veltan = temporal - repmat(sum(temporal.*geom.normal,2),[1 3]).*geom.normal;
       veladapt0 = meshadapt2(geom,parms);
       % veladapt0 = veladapt0 - repmat(sum(veladapt0.*geom.normal,2),[1 3]).*geom.normal; 
       velnormal0 = repmat(sum(velnode.*geom.normal,2),[1 3]).*geom.normal;
       if parms.maran.rkmaran ~= 0
       ajimat = surfactants(geom,velnode,veladapt0,pe,surfopt);
       end
       f4 = (velnormal0 + veladapt0 + veltan);
    elseif surfopt.opt == 3
       % hidrodinamica
       if parms.maran.rkmaran ~= 0
       ajimat = surfactants(geom,velnode,zeros(numnodes,3),pe,surfopt);
       end
       f4 = velnode0;
    elseif surfopt.opt == 4
       % hidrodinamica + adaptacion      
       veladapt0 = meshadapt2(geom,parms);
       if parms.maran.rkmaran ~= 0
       ajimat = surfactants(geom,velnode,veladapt0,pe,surfopt);
       end
       f4 = (velnode + veladapt0);
    elseif surfopt.opt == 5
           % passive (zinchenco et al. 1997)
           velnormal = repmat(sum(velnode.*geom.normal,2),[1 3]).*geom.normal;
           veladapt = meshadaptgrad(geom,velnormal,veladapt);
       if parms.maran.rkmaran ~= 0
         ajimat = surfactants(geom,velnode,veladapt,pe,surfopt);
       end
       f4 = (velnormal + veladapt);
    end
   
    if parms.maran.rkmaran ~= 0
        % evolucion del sulfactante en t + dt/2
        % propuesto
        %     gammaori = flds.gamma;
        flds.gamma = thetamethod(ajimat,deltat/2,theta,flds.gamma);

        % escale la concentracion
        gammatotnew = inttrapecioa(geom.dsi,flds.gamma)./geom.s;
        gammasc = geom.gammatotori/gammatotnew;
        flds.gamma = flds.gamma.*gammasc;
        geom.gammatot = inttrapecioa(geom.dsi,flds.gamma)./geom.s;
        geom.gammatot
    end

       geom.nodes = nodes0 + deltat*(f1+2*f2+2*f3+f4)/6;

   else
       
       % Ahora podemos usar el m?todo Predictor-Corrector de
       % Adams-Bashforth-Moulton
       
       % Paso predictor
       
       % Calculo de la velocidad en el punto actual
       
       [velnode1,geom] = stokessurf(geom,parms,flds);
       % invoque la adaptacion de la malla
       % veladapt0 = meshadapt(geom,adaptparms,velnode0);

       % factor k1 de runge kutta y calculo de las matrices de concentracion
    if surfopt.opt == 1
       % normal
       velnormal0 = repmat(sum(velnode1.*geom.normal,2),[1 3]).*geom.normal;
       if parms.maran.rkmaran ~= 0
       ajimat = surfactants(geom,velnode1,zeros(numnodes,3),pe,surfopt);
       end
       f1 = velnormal0;
    elseif surfopt.opt == 2
       % normal + adaptacion
       temporal = repmat(geom.velcentroid,[numnodes 1]);
       veltan = temporal - repmat(sum(temporal.*geom.normal,2),[1 3]).*geom.normal;
       veladapt0 = meshadapt2(geom,parms);
       % veladapt0 = veladapt0 - repmat(sum(veladapt0.*geom.normal,2),[1 3]).*geom.normal; 
       velnormal0 = repmat(sum(velnode1.*geom.normal,2),[1 3]).*geom.normal;
       if parms.maran.rkmaran ~= 0
       ajimat = surfactants(geom,velnode1,veladapt0,pe,surfopt);
       end
       f1 = (velnormal0 + veladapt0 + veltan);
    elseif surfopt.opt == 3
       % hidrodinamica
       if parms.maran.rkmaran ~= 0
       ajimat = surfactants(geom,velnode1,zeros(numnodes,3),pe,surfopt);
       end
       f1 = velnode1;
    elseif surfopt.opt == 4
       % hidrodinamica + adaptacion      
       veladapt0 = meshadapt2(geom,parms);
       if parms.maran.rkmaran ~= 0
       ajimat = surfactants(geom,velnode1,veladapt0,pe,surfopt);
       end
       f1 = (velnode1 + veladapt0);
    elseif surfopt.opt == 5
           % passive (zinchenco et al. 1997)
           velnormal = repmat(sum(velnode1.*geom.normal,2),[1 3]).*geom.normal;
           veladapt = meshadaptgrad(geom,velnormal,veladapt);
       if parms.maran.rkmaran ~= 0
         ajimat = surfactants(geom,velnode1,veladapt,pe,surfopt);
       end
       f1 = (velnormal + veladapt);
    end
   
    if parms.maran.rkmaran ~= 0
        % evolucion del sulfactante en t + dt/2
        % propuesto
        %     gammaori = flds.gamma;
        flds.gamma = thetamethod(ajimat,deltat/2,theta,flds.gamma);

        % escale la concentracion
        gammatotnew = inttrapecioa(geom.dsi,flds.gamma)./geom.s;
        gammasc = geom.gammatotori/gammatotnew;
        flds.gamma = flds.gamma.*gammasc;
        geom.gammatot = inttrapecioa(geom.dsi,flds.gamma)./geom.s;
        geom.gammatot
    end

       abm{4} = f1;
       nodes0 = geom.nodes;
       geom.nodes = geom.nodes + deltat*(-9*abm{1}+37*abm{2}-59*abm{3}+55*abm{4})/24;
       
       % Paso corrector
       
       % Calculo de la velocidad en el punto siguiente con la prediccion
       
       [velnode,geom] = stokessurf(geom,parms,flds);
       % invoque la adaptacion de la malla
       % veladapt0 = meshadapt(geom,adaptparms,velnode0);

       % factor k1 de runge kutta y calculo de las matrices de concentracion
    if surfopt.opt == 1
       % normal
       velnormal0 = repmat(sum(velnode.*geom.normal,2),[1 3]).*geom.normal;
       if parms.maran.rkmaran ~= 0
       ajimat = surfactants(geom,velnode,zeros(numnodes,3),pe,surfopt);
       end
       f2 = velnormal0;
    elseif surfopt.opt == 2
       % normal + adaptacion
       temporal = repmat(geom.velcentroid,[numnodes 1]);
       veltan = temporal - repmat(sum(temporal.*geom.normal,2),[1 3]).*geom.normal;
       veladapt0 = meshadapt2(geom,parms);
       % veladapt0 = veladapt0 - repmat(sum(veladapt0.*geom.normal,2),[1 3]).*geom.normal; 
       velnormal0 = repmat(sum(velnode.*geom.normal,2),[1 3]).*geom.normal;
       if parms.maran.rkmaran ~= 0
       ajimat = surfactants(geom,velnode,veladapt0,pe,surfopt);
       end
       f2 = (velnormal0 + veladapt0 + veltan);
    elseif surfopt.opt == 3
       % hidrodinamica
       if parms.maran.rkmaran ~= 0
       ajimat = surfactants(geom,velnode,zeros(numnodes,3),pe,surfopt);
       end
       f2 = velnode0;
    elseif surfopt.opt == 4
       % hidrodinamica + adaptacion      
       veladapt0 = meshadapt2(geom,parms);
       if parms.maran.rkmaran ~= 0
       ajimat = surfactants(geom,velnode,veladapt0,pe,surfopt);
       end
       f2 = (velnode + veladapt0);
    elseif surfopt.opt == 5
           % passive (zinchenco et al. 1997)
           velnormal = repmat(sum(velnode.*geom.normal,2),[1 3]).*geom.normal;
           veladapt = meshadaptgrad(geom,velnormal,veladapt);
       if parms.maran.rkmaran ~= 0
         ajimat = surfactants(geom,velnode,veladapt,pe,surfopt);
       end
       f1 = (velnormal + veladapt);
    end
   
    if parms.maran.rkmaran ~= 0
        % evolucion del sulfactante en t + dt/2
        % propuesto
        %     gammaori = flds.gamma;
        flds.gamma = thetamethod(ajimat,deltat/2,theta,flds.gamma);

        % escale la concentracion
        gammatotnew = inttrapecioa(geom.dsi,flds.gamma)./geom.s;
        gammasc = geom.gammatotori/gammatotnew;
        flds.gamma = flds.gamma.*gammasc;
        geom.gammatot = inttrapecioa(geom.dsi,flds.gamma)./geom.s;
        geom.gammatot
    end

       abm{5} = f2;
       
       geom.nodes = nodes0 + deltat*(abm{2}-5*abm{3}+19*abm{4}+9*abm{5})/24;
       
       % Actualizacion nodos
       
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
    paropt.tipo = 'extended';
    [geom.curv,geom.normal,geom.Kg] = curvparaboloid(geom,paropt);
    
    errorvol = abs((geom.vol - geom.volini)./geom.volini);

    for r=1:geom.numdrops
       if errorvol(r) > errorvoltol
           % invoque escalaje
           disp('Realizando Escalaje')
           geom = scaling(geom,optesc,r);
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
    
    if parms.maran.rkmaran ~= 0
   % grafique la geometria
       figure(1);
   %     grafscfld(geom,flds.gamma);
       grafscfld(geom,geom.rsigmavar);
       axis equal; view(90,0); xlabel('x1'); ylabel('x2'); zlabel('x3'); colorbar;
       getframe; title('\sigma');
       figure(2);
   %     grafscfld(geom,flds.gamma);
       grafscfld(geom,flds.gamma);
       axis equal; view(90,0); xlabel('x1'); ylabel('x2'); zlabel('x3'); colorbar;
       getframe; title('\gamma');
    else
      figure(1);
      grafscfld(geom,geom.curv);
      axis equal; view(90,0); xlabel('x1'); ylabel('x2'); zlabel('x3'); colorbar;
      getframe; title('curv');
    end
% % grafique la velocidad del centroide
%     figure(2); plot(geom.tiempo,geom.velcentroid(1,2),'*');hold on; title('velcentroid');
% grafique la velocidad normal
    figure(3); plot(geom.tiempo,velcont,'*');hold on; title('velnormal');
% grafique error de volumen
%    figure(4); plot(geom.tiempo,errorvol,'*');hold on;title('errorvol');
% % grafique la velocidad del centroide en x3
%     figure(5); plot(geom.tiempo,geom.velcentroid(1,3),'*');hold on; title('velcentroid in x3');
% % grafique la velocidad del centroide en x3
%     figure(6); plot(geom.tiempo,max(flds.gamma),'*');hold on; title('max gamma in x3');    
% % monitoreo de la concentracion
%     figure(7); plot(geom.tiempo,geom.gammatot,'*');hold on; title('gammatot');
%     
%     figure(8); plot(geom.tiempo,min(flds.gamma),'*');hold on; title('min gamma');  
%
if adim == 2
   
   % grafique la velocidad del centroide en x3
       figure(4); plot(geom.tiempo,geom.velcentroid(:,3),'*');hold on; title('vc in x3');
   % grafique la posicion del centroide en x3
       figure(5); plot(geom.tiempo,geom.xc(:,3),'*');hold on; title('xc x3');   
end    

% guarde resultados
    if counter == outputfreq
        itsaved = itsaved + 1;
        counter = 0;
        nombrearchivo = [direcciondestino num2str(itsaved), '.mat'];
        save(nombrearchivo,'geom','velnode','parms','adim','flds');
    end
    disp(carpetadestino)
    disp(['deltat: ', num2str(geom.deltat)])
    toc
end
