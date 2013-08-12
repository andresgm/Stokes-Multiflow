% clear;clc;
sbar = systembar();
% Carpeta y nombre de archivo de origen
% nombreorigenv =  {'it' 'it' 'it' 'it' 'it' 'it' 'it' 'it' 'it' 'it' 'it'
% 'it'};

carpetaorigenv = 'sed_vesicula_g01_fig4';
nombreorigenv =  'it';

% COLOQUE AQUI LA ITERACION MAXIMA QUE HAY DE CADA CARPETA DE ORIGEN
% itmaxv = [81 90 135 135 135 135 150 150 145 155 150 160];

itmaxv = 3700;

% COLOQUE AQUI LA ITERACION MINIMA (1) DE CADA CARPETA DE ORIGEN
% itminv = ones(1,size(itmaxv,2));

itminv = 1;
%itminv = ones(1,size(itmaxv,2));

% ESCRIBA AQUI EL INTERVALO DE CADA CUANTO QUIERE POPROCESAR (PILAS ESTE
% VALOR DEPENDE DE CUANTAS ITERACIONES HAY DISPONIBLE SEN LA CARPETA
% intervalv = [2 2 2 2 2 2 2 2 2 2 2 2];

intervalv = 100;

% OJO LA CANTIDAD DE ELEMENTOS DE CARPETAORIGENV, ... HASTA INTERVAL V DEBE
% SER EL MISMO.... EL RESTO ES CORRER Y YA... LE GENERA LAS IMAGENES EN
% *.FIG, *.PDF, Y *PNG DE TODO.

% OJO depende de la simulacion (sedimentacion o cortante)
ejes = [-1.5 1.5 -1.5 1.5 0 6];

nombreorigen = nombreorigenv;
carpetaorigen = carpetaorigenv;
itmin = itminv;
interval = intervalv;
itmax = itmaxv;

contador = 0;
clear pelicula xvert sigmav excesarea minxvert ...
    velcentx3 velcentx2 velcentm tiempov
%pelicula = 0;

nameydir = [cd  sbar '..' sbar 'data' sbar carpetaorigen sbar];

for k = itmin:interval:itmax
    contador = contador + 1;
    direccion = ...
        [cd  sbar '..' sbar 'data' sbar carpetaorigen sbar nombreorigen num2str(k) '.mat'];
    load(direccion);

% grafique la geometria
    figure(1);
    grafscfld(geom,geom.curv); 
    axis equal; view(90,0); 
    xlabel('x1'); ylabel('x2'); zlabel('x3'); colorbar;
    axis(ejes);
    saveas(1,[nameydir 'montage' num2str(contador)],'png')
end
