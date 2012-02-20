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

carpetaorigenv = {'pruebacortante_95_mu0_ca1'};
nombreorigenv =  {'it'};
itmaxv = [85];
itminv = ones(1,size(itmaxv,2));
intervalv = [1];

%carpetaorigenv = {'it4ele'};
%nombreorigenv =  {'it'};
%itmaxv = [69];
%itminv = [1];
%intervalv = [1];

% OJO LA CANTIDAD DE ELEMENTOS DE CARPETAORIGENV, ... HASTA INTERVAL V DEBE
% SER EL MISMO.... EL RESTO ES CORRER Y YA... LE GENERA LAS IMAGENES EN
% *.FIG, *.PDF, Y *PNG DE TODO.

% OJO dependiendo de la simulacion (sedimentacion o cortante) usted coloca los ejes
ejes = [-1.8 1.8 -1.8 1.8 -1.8 1.8];

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
        
        if strcmp(parms.flow,'inf') == 1
            velcont(contador) = max(abs(sum(velnode.*geom.normal,2)));
        end
        
        % vector de tiempo
        tiempov(contador) = geom.tiempo;

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
        
        
        figure(2); plot(tiempov,excesarea); title('excess area vs dimensionless time')
        xlabel('Dimensionless time'); ylabel('excess area');
        saveas(2,[nameydir 'excessarea'],'fig')
        saveas(2,[nameydir 'excessarea'],'eps')
        saveas(2,[nameydir 'excessarea'],'pdf')
        
         
        figure(3); plot(tiempov,def); title('Deformation vs dimensionless time');
        xlabel('Dimensionless time'); ylabel('DF = (L - B)/(L + B)');
        saveas(3,[nameydir 'def'],'fig')
        saveas(3,[nameydir 'def'],'eps')
        saveas(3,[nameydir 'def'],'pdf')    
        
        if strcmp(parms.flow,'inf') == 1
            % grafique velocidad normal
            figure(4); plot(tiempov,velcont); title('Normal Velocity vs dimensionless time');
            xlabel('Dimensionless time'); ylabel('Normal Velocity');
            saveas(4,[nameydir 'velnorm'],'fig')
            saveas(4,[nameydir 'velnorm'],'eps')
            saveas(4,[nameydir 'velnorm'],'pdf')                           
        end
        
        if thetacal == 1
            figure(5); plot(tiempov,theta); title('\theta vs dimensionless time');
            xlabel('Dimensionless time'); ylabel('\theta');
            saveas(5,[nameydir 'theta'],'fig')
            saveas(5,[nameydir 'theta'],'eps')
            saveas(5,[nameydir 'theta'],'pdf')
        
            figure(6); plot(tiempov,theta45); title('45 - \theta vs dimensionless time');
            xlabel('Dimensionless time'); ylabel('45 - \theta');
            saveas(6,[nameydir 'theta45'],'fig')
            saveas(6,[nameydir 'theta45'],'eps')
            saveas(6,[nameydir 'theta45'],'pdf')
        end
        
end
    
% valor de sigma
disp('deformation');def(end)
% exceso de area
disp('excesarea');excesarea(end)

