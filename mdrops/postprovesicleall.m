% clear;clc;
clear sigmav excesarea minxvert velcentx3 sedrate;
sbar = systembar();
% Carpeta y nombre de archivo de origen
% nombreorigenv =  {'it' 'it' 'it' 'it' 'it' 'it' 'it' 'it' 'it' 'it' 'it'
% 'it'};

carpetaorigenv = {'rbc_bagchi_relax_ref3'};
nombreorigenv =  {'it'};

% COLOQUE AQUI LA ITERACION MAXIMA QUE HAY DE CADA CARPETA DE ORIGEN
% itmaxv = [81 90 135 135 135 135 150 150 145 155 150 160];

itmaxv = 137;

% COLOQUE AQUI LA ITERACION MINIMA (1) DE CADA CARPETA DE ORIGEN
% itminv = ones(1,size(itmaxv,2));

itminv = 1;
%itminv = ones(1,size(itmaxv,2));

% ESCRIBA AQUI EL INTERVALO DE CADA CUANTO QUIERE POPROCESAR (PILAS ESTE
% VALOR DEPENDE DE CUANTAS ITERACIONES HAY DISPONIBLE SEN LA CARPETA
% intervalv = [2 2 2 2 2 2 2 2 2 2 2 2];

intervalv = [2];

% OJO LA CANTIDAD DE ELEMENTOS DE CARPETAORIGENV, ... HASTA INTERVAL V DEBE
% SER EL MISMO.... EL RESTO ES CORRER Y YA... LE GENERA LAS IMAGENES EN
% *.FIG, *.PDF, Y *PNG DE TODO.

% OJO depende de la simulacion (sedimentacion o cortante)
ejes = [-1.5 1.5 -1.5 1.5 0 9];
ejesequil = [-1.5 1.5 -1.5 1.5 0 2.5];

% Quiere que se guarde la pelicula de la sedimentacion? Si = 1, No = 0.
peliculaopt = 0;

pelicula = 0;

% Quiere la secuencia o solo los datos finales en equilibrio?
% TODO

raiz = cd;
for i = 1:size(itmaxv,2)
    nombreorigen = nombreorigenv{i};
    carpetaorigen = carpetaorigenv{i};
    itmin = itminv(i);
    interval = intervalv(i);
    itmax = itmaxv(i);
    
    contador = 0;
    clear pelicula xvert sigmav excesarea minxvert ...
        velcentx3 velcentx2 velcentm tiempov
    %pelicula = 0;
    
    for k = itmin:interval:itmax
        contador = contador + 1;
        direccion = ...
            [cd  sbar carpetaorigen sbar nombreorigen num2str(k) '.mat'];
        load(direccion);

%         % distancia del centroide a la pared
%         xvert(contador) = geom.xc(1,3);
%         % valor de sigma
%         sigmav(contador) = parms.bending.sigma;
        % exceso de area
        excesarea(contador) = abs(geom.s - geom.areaini)/geom.areaini;
%         % posicion mas baja de la gota
%         minxvert(contador) = min(geom.nodes(:,3));
%         % calculo de las deformaciones de la vesicula
%         [inerttensor,def(contador),v] = dirprindef(geom);    
        
%         % velocidad (rapidez) del centroide de la gota
%         velcentx3(contador) = abs(geom.velcentroid(1,3));
%         velcentx2(contador) = abs(geom.velcentroid(1,2));
%         velcentm(contador) = normesp(geom.velcentroid);
        % vector de tiempo
        tiempov(contador) = geom.tiempo;
        % Vector iteracion
        iterv(contador) = k;

    % grafique la geometria
        figure(1);
        grafscfld(geom,geom.curv); 
        axis equal; view(90,0); 
        xlabel('x1'); ylabel('x2'); zlabel('x3'); colorbar;
        axis(ejes);
        if peliculaopt ==1
            pelicula(contador) = getframe;
        end
    end
    
    nameydir = [raiz sbar carpetaorigen sbar];
    
    if peliculaopt == 1
        movie2avi(pelicula,[nameydir carpetaorigen '.avi'])
    end

    figure(1);
    grafscfld(geom,geom.curv);
    axis equal; view(90,0);
    xlabel('x1'); ylabel('x2'); zlabel('x3'); colorbar;
    axis([ejesequil]);
    saveas(1,[nameydir 'steadyshape'],'pdf')

%     figure(2); plot(tiempov,xvert);
%     xlabel('Time (t)');
%     ylabel('Centroid distance from wall (x_c)');
%     saveas(2,[nameydir 'xvert'],'fig')
%     saveas(2,[nameydir 'xvert'],'pdf')

%     figure(3); plot(tiempov,sigmav);
%     xlabel('Time (t)'); ylabel('\Sigma');
%     saveas(3,[nameydir 'sigma'],'fig')
%     saveas(3,[nameydir 'sigma'],'pdf')

    figure(4); plot(tiempov,excesarea);
    xlabel('Time (t)'); ylabel('Excess area (\alpha)');
    saveas(4,[nameydir 'excessarea'],'fig')
    saveas(4,[nameydir 'excessarea'],'pdf')

%     figure(5); plot(tiempov,minxvert);
%     xlabel('Dimensionless time');
%     ylabel('lowest surface vesicle position');
%     saveas(5,[nameydir 'minxvert'],'fig')
%     saveas(5,[nameydir 'minxvert'],'pdf')

%     figure(6); plot(tiempov,velcentx3);
%     xlabel('Time (t)');
%     ylabel('Vertical centroid velocity (V_{xc})');
%     saveas(6,[nameydir 'velcentx3'],'fig')
%     saveas(6,[nameydir 'velcentx3'],'pdf')        

%     figure(7); plot(tiempov,velcentx2);
%     xlabel('Dimensionless time');
%     ylabel('Centroid Velocity in X2 direction');
%     saveas(7,[nameydir 'velcentx2'],'fig')
%     saveas(7,[nameydir 'velcentx2'],'pdf')

%     figure(8); plot(tiempov,velcentm);
%     xlabel('Dimensionless time'); ylabel('Magnitude of Centroid Velocity');
%     saveas(8,[nameydir 'velcentm'],'fig')
%     saveas(8,[nameydir 'velcentm'],'pdf')

%     figure(9); plot(minxvert,velcentx3);
%     xlabel('Vesicle distance from wall');
%     ylabel('Vertical centroid velocity (V_{xc})');
%     saveas(9,[nameydir 'velcentm'],'fig')
%     saveas(9,[nameydir 'velcentm'],'pdf')
%     
%     figure(10); plot(tiempov,minxvert);
%     xlabel('Time (t)');
%     ylabel('Vesicle distance from wall');
%     saveas(10,[nameydir 'disttime'],'fig')
%     saveas(10,[nameydir 'disttime'],'pdf')
    
%     figure(11); plot(tiempov,iterv);
%     xlabel('Time (t)');
%     ylabel('Iteration');
    
%     % distancia del centroide a la pared
%     disp(['xvert = ' num2str(xvert(end))])
%     % valor de sigma
%     disp(['sigma = ' num2str(sigmav(end))])
    % exceso de area
    disp(['excesarea = ' num2str(excesarea(end))])
    % posicion mas baja de la gota
%     disp(['minxvert = ' num2str(minxvert(end))])
%     % velocidad (rapidez) del centroide de la gota
%     disp(['velcent horizontal = ' num2str(velcentx2(end))])
end
