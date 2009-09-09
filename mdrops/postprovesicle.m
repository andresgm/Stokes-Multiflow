clear;clc
% opciones de carga de archivos
    % nombre de archivo a cargar y carpeta
nombreorigen = 'it';
carpetaorigen = 'eleinf1';
itmin = 1;
interval = 1;
itmax = 74;

contador = 0;
sizevect = length(itmin:interval:itmax);
xvert = zeros(sizevect,1);
minxvert = zeros(sizevect,1);
sigmav = zeros(sizevect,1);
excesarea = zeros(sizevect,1);
velcentx3 = zeros(sizevect,1);
velcentx2 = zeros(sizevect,1);
velcentm = zeros(sizevect,1);

tiempov = zeros(sizevect,1);

nombresorigen = ['it'];

pelicula = 0;
for i = 1:size(itmaxv,2)
    nombreorigen = 'it';
    carpetaorigen = 'eleinf1';
    itmin = 1;
    interval = 1;
    itmax = 74;
    contador = 0;
    clear pelicula;
    for k = itmin:interval:itmax
        contador = contador + 1;
        direccion = [cd  '/' carpetaorigen '/' nombreorigen num2str(k) '.mat'];
        load(direccion);

        % distancia del centroide a la pared
        xvert(contador) = geom.xc(1,3);
        % valor de sigma
        sigmav(contador) = parms.bending.sigma;
        % exceso de area
        excesarea(contador) = abs(geom.s - geom.areaini)/geom.areaini;
        % posicion mas baja de la gota
        minxvert(contador) = min(geom.nodes(:,3));
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
        axis([-2 2 -2 2 0 5]);
        pelicula(contador) = getframe; title('curv');

    end
    
figure(1);
grafscfld(geom,geom.curv);
axis equal; view(90,0); xlabel('x1'); ylabel('x2'); zlabel('x3'); colorbar;
getframe; title('steady state vesicle shape - mean curvature contours');
figure(2); plot(tiempov,xvert); title('centroid position in X3 direction vs dimensionless time');
xlabel('Dimensionless time'); ylabel('centroid position');
figure(3); plot(tiempov,sigmav); title('sigma')
xlabel('Dimensionless time'); ylabel('sigma');
figure(4); plot(tiempov,excesarea); title('excess area vs dimensionless time')
xlabel('Dimensionless time'); ylabel('excess area');
figure(5); plot(tiempov,minxvert); title('lowest surface vesicle position vs dimensionless time')
xlabel('Dimensionless time'); ylabel('lowest surface vesicle position');
figure(6); plot(tiempov,velcentx3); title('Centroid Velocity in X3 direction vs dimensionless time');
xlabel('Dimensionless time'); ylabel('Centroid Velocity in X3 direction');
figure(7); plot(tiempov,velcentx2); title('Centroid Velocity in X2 direction vs dimensionless time');
xlabel('Dimensionless time'); ylabel('Centroid Velocity in X2 direction');
figure(8); plot(tiempov,velcentm); title('Magnitude of Centroid Velocity vs dimensionless time');
xlabel('Dimensionless time'); ylabel('Magnitude of Centroid Velocity');


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
