clear; clc; sbar = systembar(); raiz = cd;
% RHEOLOGY
% alpha=1/(4*pi/3);

% indique los nombres de las carpeta de origen de la soluciones Indique las
% soluciones en oreden ascendente desde el mas bajo al mas alto

carpetaorigen = 'rbc_lambda20_ca10';

% indique las iteraciones totales de la simulacion y el intervalo para los
% calculos
itmax = 999;
interval = 50;

% nombre del archivo de origen
nombreorigen = 'it';

% % Ejes [xmin xmax ymin ymax] Stress Tensor
% ejes1 = [0 10 -0.5 3.0];
% % Ejes [xmin xmax ymin ymax] Visc Particle

contador = 0;
y1 = []; y2 = []; y3 = []; y4 = []; time = [];

for k = 1:interval:itmax
    
    contador = contador + 1;
    % cargue el archivo de solucion
    direccion = ...
        [cd  sbar '..' sbar 'data' sbar carpetaorigen sbar ...
        nombreorigen num2str(k) '.mat'];
    load(direccion);
   
    % calcule el delta de fuerza en la interfase por curvatura
    rdeltaf = geom.rdeltafnorm+geom.rdeltafmaran;
    
    strtensornode = strtensorreologia(rdeltaf,geom.nodes,...
        velnode,geom.normal,parms.lamda);
    rintegral = inttrapeciomata(geom.dsi,strtensornode)/geom.vol;
    
    time(contador) = geom.tiempo;
    % S11 - S22 (S22 - S33 este caso)  
    y1(contador) = rintegral(2,2)-rintegral(3,3);
    % S22 - S33 (S33 - S11 este caso)
    y2(contador) = rintegral(3,3)-rintegral(1,1);
    % S12 (S23 este caso)
    y3(contador) = rintegral(2,3);
    % Sii/3 - 2
    y4(contador) = trace(rintegral)/3 - 2;

end

figure(1);
Name = ['Stress Tensor for \lambda=' num2str(parms.lamda)...
    ' and Ca=' num2str(parms.ca)];
plot(time,y1,'--k',time,y2,'-.k',time,y3,'-k',time,y4,':k','LineWidth',1)
xlabel('Time','fontsize',14)
ylabel('\tau_{ij}','fontsize',14)
title(Name,'fontsize',14);
leyenda=legend('\tau_{11}-\tau_{22}','\tau_{22}-\tau_{33}',...
    '\tau_{12}','(\tau_{ii}/{3})-2',2);
set(leyenda,'fontsize',14)

nameydir = [raiz sbar];
% axis([ejes1]);
saveas(1,[nameydir 'strtenlamda' num2str(parms.lamda*10)],'fig')
saveas(1,[nameydir 'strtenlamda' num2str(parms.lamda*10)],'eps')
saveas(1,[nameydir 'strtenlamda' num2str(parms.lamda*10)],'pdf')
        
% nondimensionless viscosity pg 111 tesis kennedy
y7 = y3(1:end)./parms.ca;

figure(2);
Name = ['Viscosity \lambda=' num2str(parms.lamda) ' and Ca=' num2str(parms.ca)];
plot(time(1:end),y7,'-k','LineWidth',0.5)
xlabel('Time','fontsize',14)
ylabel('\mu*=\tau_{12}/ca','fontsize',14)
set(gca,'fontsize',14)
title(Name,'fontsize',14);
grid on;

nameydir = [raiz sbar];
ejes2 = [0 geom.tiempo 0 2.5];
axis([ejes2]);
saveas(2,[nameydir 'visclamda' num2str(parms.lamda*10)],'fig')
saveas(2,[nameydir 'visclamda' num2str(parms.lamda*10)],'eps')
saveas(2,[nameydir 'visclamda' num2str(parms.lamda*10)],'pdf')