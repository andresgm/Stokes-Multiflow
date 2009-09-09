clear;clc
sbar = systembar();
% COLOQUE AQUI LA ITERACION MAXIMA QUE HAY DE CADA CARPETA DE ORIGEN

% COLOQUE AQUI LA ITERACION MINIMA (1) DE CADA CARPETA DE ORIGEN

% ESCRIBA AQUI EL INTERVALO DE CADA CUANTO QUIERE POPROCESAR (PILAS ESTE
% VALOR DEPENDE DE CUANTAS ITERACIONES HAY DISPONIBLE SEN LA CARPETA

% nombreorigenv =  {'it' 'it' 'it' 'it' 'it' 'it' 'it' 'it' 'it' 'it' 'it' 'it'};
% itmaxv = [81 90 135 135 135 135 150 150 145 155 150 160];
% itminv = ones(1,size(itmaxv,2));
% intervalv = [2 2 2 2 2 2 2 2 2 2 2 2];

% carpetaorigenv = {'shearinf2' 'shearinf3' 'shearinf4' 'shearinf5' 'shearinf6' };
% nombreorigenv =  {'it' 'it' 'it' 'it' 'it'};
% itmaxv = [138 354 358 166 153];
% itminv = ones(1,size(itmaxv,2));
% intervalv = [1 2 2 1 1];

carpetaorigenv = {'gotaLambda3.6Ca1.5'};
nombreorigenv =  {'it'};
itmaxv = [600];
itminv = ones(1,size(itmaxv,2));
intervalv = [5];

%carpetaorigenv = {'it4ele'};
%nombreorigenv =  {'it'};
%itmaxv = [69];
%itminv = [1];
%intervalv = [1];

% OJO LA CANTIDAD DE ELEMENTOS DE CARPETAORIGENV, ... HASTA INTERVAL V DEBE
% SER EL MISMO.... EL RESTO ES CORRER Y YA... LE GENERA LAS IMAGENES EN
% *.FIG, *.PDF, Y *PNG DE TODO.

% OJO dependiendo de la simulacion (sedimentacion o cortante) usted coloca los ejes
ejes = [-3 3 -3 3 -3 3];

% voy a ir a almorzar mientras posprocese lo que hace falta

% calcule el angulo theta? si:1 no:0
thetacal = 1;

pelicula = 0;
raiz = cd;
for i = 1:size(itmaxv,2)
    nombreorigen = nombreorigenv{i};
    carpetaorigen = carpetaorigenv{i};
    itmin = itminv(i);
    interval = intervalv(i);
    itmax = itmaxv(i);
    
    contador = 0;
    clear pelicula xvert sigmav excesarea minxvert velcentx3 velcentx2 velcentm tiempov
    
    for k = itmin:interval:itmax
        contador = contador + 1;
        direccion = [cd  sbar carpetaorigen sbar nombreorigen num2str(k) '.mat'];
        load(direccion);

        % distancia del centroide a la pared
        xvert(contador) = geom.xc(1,3);
        % exceso de area
        excesarea(contador) = abs(geom.s - geom.areaini)/geom.areaini;
        % posicion mas baja de la gota
        minxvert(contador) = min(geom.nodes(:,3));
        % calculo de las deformaciones de la gota
        if thetacal == 1
            [inerttensor,def(contador),v,theta(contador)] = dirprindef(geom,1);    
            theta45(contador) = 45 - theta(contador);
        else
            [inerttensor,def(contador),v,theta] = dirprindef(geom);    
        end
        
        % velocidad (rapidez) del centroide de la gota
        velcentx3(contador) = abs(geom.velcentroid(1,3));
        velcentx2(contador) = abs(geom.velcentroid(1,2));
        velcentm(contador) = normesp(geom.velcentroid);
        if strcmp(parms.flow,'inf') == 1
            velcont(contador) = max(abs(sum(velnode.*geom.normal,2)));
        end
        % vector de tiempo
        tiempov(contador) = geom.tiempo;

        geom.deltat

    % grafique la geometria
        figure(1);
        grafscfld(geom,geom.curv); 
        axis equal; view(90,0); xlabel('x1'); ylabel('x2'); zlabel('x3'); colorbar;
        axis(ejes);
        pelicula(contador) = getframe; title('curv');
        
    end
        nameydir = [raiz sbar carpetaorigen sbar];
        movie2avi(pelicula,[nameydir carpetaorigen '.avi'])
        figure(1);
        grafscfld(geom,geom.curv);
        axis equal; view(90,0); xlabel('x1'); ylabel('x2'); zlabel('x3'); colorbar;
        title('steady state vesicle shape - mean curvature contours');
        axis([ejes]);
        saveas(1,[nameydir 'shape'],'fig')
        saveas(1,[nameydir 'shape'],'eps')
        saveas(1,[nameydir 'shape'],'pdf')
        
        figure(2); plot(tiempov,xvert); title('centroid position in X3 direction vs dimensionless time');
        xlabel('Dimensionless time'); ylabel('centroid position');
        saveas(2,[nameydir 'xvert'],'fig')
        saveas(2,[nameydir 'xvert'],'eps')
        saveas(2,[nameydir 'xvert'],'pdf')
        
        figure(4); plot(tiempov,excesarea); title('excess area vs dimensionless time')
        xlabel('Dimensionless time'); ylabel('excess area');
        saveas(4,[nameydir 'excessarea'],'fig')
        saveas(4,[nameydir 'excessarea'],'eps')
        saveas(4,[nameydir 'excessarea'],'pdf')
        
        figure(5); plot(tiempov,minxvert); title('lowest surface vesicle position vs dimensionless time')
        xlabel('Dimensionless time'); ylabel('lowest surface vesicle position');
        saveas(5,[nameydir 'minxvert'],'fig')
        saveas(5,[nameydir 'minxvert'],'eps')
        saveas(5,[nameydir 'minxvert'],'pdf')
        
        figure(6); plot(tiempov,velcentx3); title('Centroid Velocity in X3 direction vs dimensionless time');
        xlabel('Dimensionless time'); ylabel('Centroid Velocity in X3 direction');
        saveas(6,[nameydir 'velcentx3'],'fig')
        saveas(6,[nameydir 'velcentx3'],'eps')
        saveas(6,[nameydir 'velcentx3'],'pdf')        
        
        figure(7); plot(tiempov,velcentx2); title('Centroid Velocity in X2 direction vs dimensionless time');
        xlabel('Dimensionless time'); ylabel('Centroid Velocity in X2 direction');
        saveas(7,[nameydir 'velcentx2'],'fig')
        saveas(7,[nameydir 'velcentx2'],'eps')
        saveas(7,[nameydir 'velcentx2'],'pdf')
        
        figure(8); plot(tiempov,velcentm); title('Magnitude of Centroid Velocity vs dimensionless time');
        xlabel('Dimensionless time'); ylabel('Magnitude of Centroid Velocity');
        saveas(8,[nameydir 'velcentm'],'fig')
        saveas(8,[nameydir 'velcentm'],'eps')
        saveas(8,[nameydir 'velcentm'],'pdf')       
        
        figure(9); plot(tiempov,def); title('Deformation vs dimensionless time');
        xlabel('Dimensionless time'); ylabel('DF = (L - B)/(L + B)');
        saveas(9,[nameydir 'def'],'fig')
        saveas(9,[nameydir 'def'],'eps')
        saveas(9,[nameydir 'def'],'pdf')    
        
        if strcmp(parms.flow,'inf') == 1
            % grafique velocidad normal
            figure(10); plot(tiempov,velcont); title('Normal Velocity vs dimensionless time');
            xlabel('Dimensionless time'); ylabel('Normal Velocity');
            saveas(10,[nameydir 'velnorm'],'fig')
            saveas(10,[nameydir 'velnorm'],'eps')
            saveas(10,[nameydir 'velnorm'],'pdf')                           
        end
        
        if thetacal == 1
            figure(11); plot(tiempov,theta); title('\theta vs dimensionless time');
            xlabel('Dimensionless time'); ylabel('\theta');
            saveas(11,[nameydir 'theta'],'fig')
            saveas(11,[nameydir 'theta'],'eps')
            saveas(11,[nameydir 'theta'],'pdf')
        
            figure(12); plot(tiempov,theta45); title('45 - \theta vs dimensionless time');
            xlabel('Dimensionless time'); ylabel('45 - \theta');
            saveas(12,[nameydir 'theta45'],'fig')
            saveas(12,[nameydir 'theta45'],'eps')
            saveas(12,[nameydir 'theta45'],'pdf')
        end
        
end
    
% distancia del centroide a la pared
disp('xvert');xvert(end)
% valor de sigma
disp('deformation');def(end)
% exceso de area
disp('excesarea');excesarea(end)
% posicion mas baja de la gota
disp('minxvert');minxvert(end)
% velocidad (rapidez) del centroide de la gota
disp('velcent horizontal');velcentx2(end)
