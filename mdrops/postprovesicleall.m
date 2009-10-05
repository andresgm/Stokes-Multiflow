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

carpetaorigenv = {'sedimentacion_g0_1'};
nombreorigenv =  {'it'};
itmaxv = [8000];
itminv = ones(1,size(itmaxv,2));
intervalv = [100];

%carpetaorigenv = {'eleinf1'};
%nombreorigenv =  {'it' };
%itmaxv = [76];
%itminv = ones(1,size(itmaxv,2));
%intervalv = [1];

%carpetaorigenv = {'it4ele'};
%nombreorigenv =  {'it'};
%itmaxv = [69];
%itminv = [1];
%intervalv = [1];

% OJO LA CANTIDAD DE ELEMENTOS DE CARPETAORIGENV, ... HASTA INTERVAL V DEBE
% SER EL MISMO.... EL RESTO ES CORRER Y YA... LE GENERA LAS IMAGENES EN
% *.FIG, *.PDF, Y *PNG DE TODO.

% OJO dependiendo de la simulacion (sedimentacion o cortante) usted coloca los ejes
ejes = [-1.5 1.5 -1.5 1.5 0 11];

% voy a ir a almorzar mientras posprocese lo que hace falta


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
        % valor de sigma
        sigmav(contador) = parms.bending.sigma;
        % exceso de area
        excesarea(contador) = abs(geom.s - geom.areaini)/geom.areaini;
        % posicion mas baja de la gota
        minxvert(contador) = min(geom.nodes(:,3));
        % calculo de las deformaciones de la vesicula
        [inerttensor,def(contador),v] = dirprindef(geom);    
        
        % velocidad (rapidez) del centroide de la gota
        velcentx3(contador) = abs(geom.velcentroid(1,3));
        velcentx2(contador) = abs(geom.velcentroid(1,2));
        velcentm(contador) = normesp(geom.velcentroid);
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
        saveas(1,[nameydir 'shape'],'png')
        saveas(1,[nameydir 'shape'],'pdf')
        
        figure(2); plot(tiempov,xvert); title('centroid position in X3 direction vs dimensionless time');
        xlabel('Dimensionless time'); ylabel('centroid position');
        saveas(2,[nameydir 'xvert'],'fig')
        saveas(2,[nameydir 'xvert'],'png')
        saveas(2,[nameydir 'xvert'],'pdf')
        
        figure(3); plot(tiempov,sigmav); title('sigma')
        xlabel('Dimensionless time'); ylabel('sigma');
        saveas(3,[nameydir 'sigma'],'fig')
        saveas(3,[nameydir 'sigma'],'png')
        saveas(3,[nameydir 'sigma'],'pdf')
        
        figure(4); plot(tiempov,excesarea); title('excess area vs dimensionless time')
        xlabel('Dimensionless time'); ylabel('excess area');
        saveas(4,[nameydir 'excessarea'],'fig')
        saveas(4,[nameydir 'excessarea'],'png')
        saveas(4,[nameydir 'excessarea'],'pdf')
        
        figure(5); plot(tiempov,minxvert); title('lowest surface vesicle position vs dimensionless time')
        xlabel('Dimensionless time'); ylabel('lowest surface vesicle position');
        saveas(5,[nameydir 'minxvert'],'fig')
        saveas(5,[nameydir 'minxvert'],'png')
        saveas(5,[nameydir 'minxvert'],'pdf')
        
        figure(6); plot(tiempov,velcentx3); title('Centroid Velocity in X3 direction vs dimensionless time');
        xlabel('Dimensionless time'); ylabel('Centroid Velocity in X3 direction');
        saveas(6,[nameydir 'velcentx3'],'fig')
        saveas(6,[nameydir 'velcentx3'],'png')
        saveas(6,[nameydir 'velcentx3'],'pdf')        
        
        figure(7); plot(tiempov,velcentx2); title('Centroid Velocity in X2 direction vs dimensionless time');
        xlabel('Dimensionless time'); ylabel('Centroid Velocity in X2 direction');
        saveas(7,[nameydir 'velcentx2'],'fig')
        saveas(7,[nameydir 'velcentx2'],'png')
        saveas(7,[nameydir 'velcentx2'],'pdf')
        
        figure(8); plot(tiempov,velcentm); title('Magnitude of Centroid Velocity vs dimensionless time');
        xlabel('Dimensionless time'); ylabel('Magnitude of Centroid Velocity');
        saveas(8,[nameydir 'velcentm'],'fig')
        saveas(8,[nameydir 'velcentm'],'png')
        saveas(8,[nameydir 'velcentm'],'pdf')        
end
    
% distancia del centroide a la pared
disp('xvert');xvert(end)
% valor de sigma
disp('sigma');sigmav(end)
% exceso de area
disp('excesarea');excesarea(end)
% posicion mas baja de la gota
disp('minxvert');minxvert(end)
% velocidad (rapidez) del centroide de la gota
disp('velcent horizontal');velcentx2(end)
