% CALCULO DEL FLUJO DE STOKES PARA UNA GOTA
% IMPLEMENTADO FLUJO INFINITO, SEMIINFINITO
% IMPLEMENTADO SINGLE Y DOUBLE LAYER
clear;clc; %close all
%% opciones de carga de archivos
    % nombre de archivo a cargar y carpeta
nombreorigen = 'it';
carpetaorigen = 'drops_la_1_ca_0.04_x_0_ext';
iteracion = 10;

% nombre de archivo a guardar y carpeta
nombredestino = 'it';
carpetadestino = 'drops_la_1_ca_0.04_x_0_ext_testmeshadaptcurv';
% simulacion nueva desde cero optsim = 0
% continue la simulacion optsim = 1
% simulacion nueva desde archivo de resultados optsim = 2
opcionsim = 2;

% Algoritmo de flujo de stokes.
ca = 0.04;
lamda = 1;

% tipo de flujo flow: 'inf' flow:'semiinf'
flow = 'inf';

%Forma del flujo 0=Cortante simple, 1=plano hiperbolico ,2=Extensional Puro
Fflow = 2;

% opcion de calculo de la curvatura 1: paraboloid fitting; 2: extended par;
% 3: basado en laplace beltrami
curvopt = 3;

% Banderas de fuerza dif 0: si. 1: no
% gravedad
kb = 0;
Bo = 1;

% Tensoactivos
kc = 0;
% maranmodel = 1(lineal) definir beta y pe(alpha)
% maranmodel = 2(logaritmico) definir x, e y pe(alpha)
maranmodel = 2;
% Parametro \Alpha = \simga_0 R_0/\mu D_s
alpha = 1000;
% constante para modelo lineal
beta = 0.2;
% Parametro del modelo logaritmico
e = 0.2;
% Concentracion
x = 0.1;

%Solubilidad
ks=0;
%Parametro de solubilidad B
B=0.01;
%parametro de profundidad k
k= -x/(x-1);
%numero de Biot (Bi)
Bi=B/ca;

% numero de gotas
geom.numdrops = 1;
% Coordenadas de los centroides de las gotas
xc =[0 0 0];

% frecuencia de guardar resultados
outputfreq = 10;

% pasos de tiempo de la simulacion
numtimesteps = 80000;
% Reduccion del paso de tiempo calculado automaticamente
% redfactor = 5000;
deltat = 0.1*ca;

% escalaje
errorvoltol = 1e-6;
optesc.maxit = 15000;
optesc.kp = 20;
optesc.deltate = 0.01;
optesc.tolerrorvol = errorvoltol;

% parametros de adaptacion
% velopt: 1 hidrodinamica velopt:2 normal velopt:3 passive (zinchenko et al.)
velopt = 3;
surfopt.opt = 5;

% parametros de tiempo para integracion de surfactantes
% theta = 0 Euler Explicito
% theta = 0.5 semi implicito
% theta = 1 full implicito
theta = 1;

%% procesamiento de parametros
% adimensionalizacion de andres gonzalez basado en la velocidad
% caracteristica de sedimentacion U_0=\Delata \rho g R_0^2/\mu(1+\lambda)
parms.flow = flow;

% adimensionalizacion del single layer
parms.rkextf = 2*ca;
parms.rksl = 2;
parms.rkdl = 2*(lamda - 1)/(lamda + 1);
parms.lamda = lamda;
parms.ca = ca;
parms.Bo = Bo;
parms.Fflow = Fflow;
parms.Bi= Bi;
parms.k= k;

% Coeficiente termino de curvatura
parms.rkcurv = 1;

if kb == 1
    % gravedad
    parms.rkgrav = Bo;
else
    parms.rkgrav = 0;
end

if kc == 1
    gammaeqovergamma0 = 1+e*log(1-x);
    pe = ca*alpha*gammaeqovergamma0;
    parms.maran.rkmaran = 1;
    if maranmodel == 1
        parms.maran.maranmodel = 'linear';
        parms.maran.beta = beta;
        parms.maran.pe = pe;
    elseif maranmodel == 2
        parms.maran.maranmodel = 'log';
        parms.maran.pe = pe;
        parms.maran.x = x;
        parms.maran.e = e;
    end
else
    parms.maran.rkmaran = 0;
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
    
    if parms.lamda ~= 0
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
    geom.normalele = geomprop.normalele;
    geom.normal = geomprop.normal;
    geom.dsi = geomprop.dsi;
    geom.ds = geomprop.ds;
    geom.s = geomprop.s;
    geom.vol = geomprop.vol;
    geom.jacmat = geomprop.jacmat;
    geom.g = geomprop.g;
    geom.volini = geom.vol;
    geom.areaini = geom.s;
    geom.xcini = centroide(geom);
    
    % Calculo inicial de la curvatura usando ajuste a superficie cuadratica
    
    paropt.tipo = 'extended';
    [geom.curv,geom.normal,geom.Kg] = curvparaboloid(geom,paropt);
    
    direcciondestino = ...
        [cd sbar '..' sbar 'data' sbar carpetadestino sbar nombredestino];
    mkdir([cd sbar '..' sbar 'data' sbar carpetadestino]);
        
    paso = 1;
    counter = 0;
    geom.tiempo = 0;
    itsaved = 0;
    
    % Campo inicial de concentracion
    flds.gamma = ones(numnodes,1);
    geom.gammatotori = inttrapecioa(geom.dsi,flds.gamma)./geom.s;

elseif opcionsim == 1
    % cargue desde resultados y continue la simulacion
    carpetadestino = carpetaorigen;
    nombredestino = nombreorigen;
    if isempty(carpetaorigen) == 1
        error('Must specify folder to read previous result.');
    else
        direccion = ...
         [cd sbar '..' sbar 'data' sbar carpetaorigen sbar nombreorigen...
         num2str(iteracion) '.mat'];
    end
    
    load(direccion);
    
    paso = iteracion + 1;
    counter = 0;
    itsaved = iteracion;

    direcciondestino = ...
        [cd sbar '..' sbar 'data' sbar carpetadestino sbar nombredestino];
    mkdir([cd sbar '..' sbar 'data' sbar carpetadestino]);
    numnodes = size(geom.nodes,1);
    numelements = size(geom.elements,1);
    
    parms.curvopt = parmstemp.curvopt;
    
elseif opcionsim == 2
    % cargue desde resultados y realice una nueva simulacion
    % cargue desde resultados y continue la simulacion
    if isempty(carpetaorigen) == 1
        error('Must specify folder to read previous result.');
    else
        direccion = ...
         [cd sbar '..' sbar 'data' sbar carpetaorigen sbar nombreorigen...
         num2str(iteracion) '.mat'];
    end
    
    load(direccion);
    
    direcciondestino = ...
        [cd sbar '..' sbar 'data' sbar carpetadestino sbar nombredestino];
    mkdir([cd sbar '..' sbar 'data' sbar carpetadestino]);
        
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
    
end
  
%% Ciclo principal
% xcant = centroide(geom);
% geom.velcentroid = [0 0 0];
% geom.xc = xcant;
defant = 100;
veladapt = zeros(numnodes,3);

for p = paso:numtimesteps
%tic
% if p == 5
% profile on
% end
% if p == 10
% profile off
% end

    counter = counter + 1;
    disp('************************')
    disp(['iteracion = ', num2str(p)])
    
% calcule la distancia minima de adaptacion y el paso de tiempo
    
% lmin = zeros(numnodes,1);
% for k = 1:numnodes
% % extraiga los nodos vecinos a un nodo en la misma gota
% nodesadj = geom.nodecon2node{k};
% lmin(k) = min(normesp(repmat(geom.nodes(k,:),[size(nodesadj,1) 1])...
% - geom.nodes(nodesadj,:)));
% end
% % Longitud minima para verificacion de paso de tiempo.
% lmint = min(lmin);
% deltat = lmint^1.5/redfactor;
% parms.lmin = lmin;
    
% Usamos esquema rk2 para la integracion numerica.
    
    %% primer paso de runge kutta f1
    
    [velnode1,geom] = stokessurf(geom,parms,flds);

    % passive (zinchenco et al. 1997)
    velnormal = repmat(sum(velnode1.*geom.normal,2),[1 3]).*geom.normal;
    veladapt = meshadaptgradcurv(geom,velnormal,veladapt,deltat/2);
    f1 = (velnormal + veladapt);
   
    if parms.maran.rkmaran ~= 0
        % Third argument is the adaptation velocity w in Bazhlekov et al
        % 2003 sense, careful. v = u + w.
         if ks==0
            
            ajimat = ...
                surfactants(geom,velnode1,veladapt-(velnode1-velnormal),pe);
        else
            ajimat = ...
                solsurf(geom,velnode1,veladapt-(velnode1-velnormal),pe,Bi,k);
        end
        % evolucion del sulfactante en t + dt/2
        % propuesto
        gammaori = flds.gamma;
        flds.gamma = thetamethod(ajimat,deltat/2,theta,flds.gamma);

        % escale la concentracion
% gammatotnew = inttrapecioa(geom.dsi,flds.gamma)./geom.s;
% gammasc = geom.gammatotori/gammatotnew;
% flds.gamma = flds.gamma.*gammasc;
% geom.gammatot = inttrapecioa(geom.dsi,flds.gamma)./geom.s;
% geom.gammatot;
    end

    nodes0 = geom.nodes;
       
    geom.nodes = geom.nodes + (1/2)*deltat*f1;
       
    %% segundo paso de runge kutta f2
    % invoque el problema de flujo de stokes
    
    % calcule el vector normal a cada nodo
    normalandgeoopt.normal = 1;
    normalandgeoopt.areas = 1;
    normalandgeoopt.vol = 1;
    geomprop = normalandgeo(geom,normalandgeoopt,1);
    geom.normalele = geomprop.normalele;
    geom.normal = geomprop.normal;
    geom.dsi = geomprop.dsi;
    geom.ds = geomprop.ds;
    geom.s = geomprop.s;
    geom.vol = geomprop.vol;
    geom.jacmat = geomprop.jacmat;
    geom.g = geomprop.g;
    paropt.tipo = 'extended';
    [geom.curv,geom.normal,geom.Kg] = curvparaboloid(geom,paropt);
    
    [velnode,geom] = stokessurf(geom,parms,flds);
    
    % passive (zinchenco et al. 1997)
    velnormal = repmat(sum(velnode.*geom.normal,2),[1 3]).*geom.normal;
    veladapt = meshadaptgradcurv(geom,velnormal,veladapt,deltat);
    f2 = (velnormal + veladapt);
   
    if parms.maran.rkmaran ~= 0
        % Third argument is the adaptation velocity w in Bazhlekov et al
        % 2003 sense, careful. v = u + w.
        if ks==0
            
            ajimat = ...
                surfactants(geom,velnode1,veladapt-(velnode1-velnormal),pe);
        else
            ajimat = ...
                solsurf(geom,velnode1,veladapt-(velnode1-velnormal),pe,Bi,k);
        end
        % evolucion del sulfactante en t
        flds.gamma = thetamethod(ajimat,deltat,theta,gammaori);

        % escale la concentracion
% gammatotnew = inttrapecioa(geom.dsi,flds.gamma)./geom.s;
% gammasc = geom.gammatotori/gammatotnew;
% flds.gamma = flds.gamma.*gammasc;
% geom.gammatot = inttrapecioa(geom.dsi,flds.gamma)./geom.s;
% geom.gammatot
    end

    geom.nodes = nodes0 + deltat*f2;
    
%% escalaje
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
    paropt.tipo = 'extended';
    [geom.curv,geom.normal,geom.Kg] = curvparaboloid(geom,paropt);
    
    errorvol = abs((geom.vol - geom.volini)./geom.volini);

    if errorvol > errorvoltol
       % invoque escalaje
       geom = scaling(geom,optesc,errorvol);
    end

% error de volumen
% errorvol = abs(geom.volini - geom.vol)/geom.volini;
% disp(['Error volumen post-escalaje: ',num2str(errorvol)]);
% velocidad normal maxima y tiempo de simulacion

    velcont = max(abs(sum(velnode.*geom.normal,2)));
    disp(['Velocidad normal maxima: ', num2str(velcont)]);
    geom.tiempo = geom.tiempo + deltat;
    geom.deltat = deltat;
    
% calculo del centroide de la gota y velocidad del centroide
% xcant = geom.xc;
% geom.xc = centroide(geom);
% geom.velcentroid = (geom.xc - xcant)./deltat;
% disp(['Posicion centroide: ', num2str([geom.xc(1),geom.xc(2),geom.xc(3)])]);
% disp(['Velocidad centroide: ', num2str(geom.velcentroid)]);
    
    if parms.maran.rkmaran ~= 0
   % grafique la geometria
       figure(1);
   % grafscfld(geom,flds.gamma);
       grafscfld(geom,normesp(geom.rdeltafcurv));
       axis equal; view(90,0); xlabel('x1'); ylabel('x2'); zlabel('x3'); colorbar;
       getframe; title('Tension');
       figure(2);
   % grafscfld(geom,flds.gamma);
       grafscfld(geom,normesp(geom.rdeltafmaran));
       axis equal; view(90,0); xlabel('x1'); ylabel('x2'); zlabel('x3'); colorbar;
       getframe; title('Marangoni');
    else
      figure(1);
      grafscfld(geom,geom.curv);
      axis equal; view(90,0); xlabel('x1'); ylabel('x2'); zlabel('x3'); colorbar;
      getframe; title('curv');
    end
    
% if p > 1
% figure(3); plot(geom.tiempo,geom.velcentroid(3),'*r');hold on;
% title('Sedimentation Rate'); getframe;
% end

% if kb == 1
% disp(['Fuerzagrav: ', num2str(geom.fuerzagrav)]);
% end

    [inerttensor,def,v,theta] = dirprindef(geom,1);
    disp(['DF: ', num2str(def)]);
    disp(['theta: ', num2str((theta)/pi)]);
    vardef = abs((defant-def)/def);
    disp(['Cambio DF: ', num2str(vardef)]);
    if vardef <= 1e-6
        disp('Convergencia en deformacion');
        break
    else
        defant = def;
    end
    
%     figure(3); plot(geom.tiempo,def,'*k');hold on;
%     title('DF'); 

% guarde resultados
    if counter == outputfreq
        itsaved = itsaved + 1;
        counter = 0;
        nombrearchivo = [direcciondestino num2str(itsaved), '.mat'];
        save(nombrearchivo,'geom','velnode','parms','flds');
    end
% disp(carpetadestino)
    disp(['deltat: ', num2str(geom.deltat)]);
    disp(['tiempo: ', num2str(geom.tiempo)]);
    
% if velcont*deltat < 1e-10
% disp('Convergencia a estado estacionario');
% itsaved = itsaved + 1;
% nombrearchivo = [direcciondestino num2str(itsaved), '.mat'];
% save(nombrearchivo,'geom','velnode','parms');
% break;
% end

% if p == 14
% profile viewer
% end
    
% toc
end