clear; clc; sbar = systembar(); raiz = cd;
% RHEOLOGY
alpha=1/(4*pi/3);

% indique los nombres de las carpeta de origen de la soluciones Indique las
% soluciones en oreden ascendente desde el mas bajo al mas alto

% % Lambda = 6.4
% carpetaorigen = {'gotaLambda6.4Ca0.05' 'gotaLambda6.4Ca0.1' 'gotaLambda6.4Ca0.25'...
%     'gotaLambda6.4Ca0.50' 'gotaLambda6.4Ca0.75' 'gotaLambda6.4Ca1.00' 'gotaLambda6.4Ca1.50'}
% 
% % indique las iteraciones de estado estable para cada solucion
% itss = [695 281 279 270 639 255 984];

% % Lambda = 3.6
% carpetaorigen = {'gotaLambda3.6Ca0.05' 'gotaLambda3.6Ca0.1' 'gotaLambda3.6Ca0.25'...
%     'gotaLambda3.6Ca0.50' 'gotaLambda3.6Ca0.75' 'gotaLambda3.6Ca1' 'gotaLambda3.6Ca1.5' ...
%     'gotaLambda3.6Ca2.00'}
% 
% % indique las iteraciones de estado estable para cada solucion
% itss = [699 1072 1068 1073 1080 1132 1042 2085];

carpetaorigen = {'gotaLambda0.2Ca0.35'}

% indique las iteraciones de estado estable para cada solucion
itss = [30];

% nombre del archivo de origen
nombreorigen = 'it';

% Ejes [xmin xmax ymin ymax] Stress Tensor
ejes1 = [0 0.45 -0.5 3.0];
% Ejes [xmin xmax ymin ymax] Visc Particle
ejes2 = [0 0.45 0 2.5];

contador = 1;
y1(1) = 0; y2(1) = 0; y3(1) = 0; y4(1) = 0; capilar(1) = 0;

for k=1: max(size(carpetaorigen))
    contador = contador + 1;
    % cargue el archivo de solucion
    direccion = [cd  sbar carpetaorigen{k} sbar nombreorigen num2str(itss(k)) '.mat'];
    load(direccion);
   
    % calcule el delta de fuerza en la interfase por curvatura
    rdeltafcurv = deltafcurv(geom.curv,parms.rkcurv);
    deltafcurvv = geom.normal.*repmat(rdeltafcurv,[1 3]);
   
    rintegral = zeros(3,3);
    for z=1:geom.numdrops
        strtensornode = strtensorreologia(deltafcurvv,geom.nodes,...
            velnode,geom.normal,parms.lamda,geom.nnodesdrop,z);
        rintegral = rintegral + alpha*inttrapeciomata(geom.dsi(geom.nnodesdrop(z,1):...
            geom.nnodesdrop(z,2),:),strtensornode);
    end
    
    capilar(contador) = parms.ca;
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
Name = ['Stress Tensor for \lambda=' num2str(parms.lamda)];
plot(capilar,y1,'-b+',capilar,y2,'-ro',capilar,y3,'-gs',capilar,y4,'-c>','LineWidth',0.1,'MarkerEdgeColor','k','MarkerSize',1)
xlabel('Ca','fontsize',12)
ylabel('\tau_{ij}','fontsize',12)
title(Name,'fontsize',12);
leyenda=legend('\tau_{11}-\tau_{22}','\tau_{22}-\tau_{33}','\tau_{12}','(\tau_{ii}/{3})-2',2);
set(leyenda,'fontsize',12)

nameydir = [raiz sbar];
axis([ejes1]);
saveas(1,[nameydir 'strtenlamda' num2str(parms.lamda*10)],'fig')
saveas(1,[nameydir 'strtenlamda' num2str(parms.lamda*10)],'eps')
saveas(1,[nameydir 'strtenlamda' num2str(parms.lamda*10)],'pdf')
        
% nondimensionless viscosity pg 111 tesis kennedy
y7 = y3(2:end)./capilar(2:end);

figure(2);
Name = ['Viscosity \lambda=' num2str(parms.lamda)];
plot(capilar(2:end),y7,'-ko','LineWidth',0.5,'MarkerEdgeColor','k','MarkerSize',2)
xlabel('Ca','fontsize',12)
ylabel('\mu*=\tau_{12}/ca','fontsize',12)
set(gca,'fontsize',12)
title(Name,'fontsize',12);
grid on;

nameydir = [raiz sbar];
axis([ejes2]);
saveas(2,[nameydir 'visclamda' num2str(parms.lamda*10)],'fig')
saveas(2,[nameydir 'visclamda' num2str(parms.lamda*10)],'eps')
saveas(2,[nameydir 'visclamda' num2str(parms.lamda*10)],'pdf')