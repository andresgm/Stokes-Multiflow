% calcula el tensor de inercial de una superficie cerrada
% calcula las coordendas de los nodos en la base ubicada en el centroide
% del cuerpo
% TODO: definir documentacion

function [inerttensor,nodesloc] = itensor(geom)

numnodes = size(geom.nodes,1);

if isfield(geom,'normal') ~= 1
    % calcule el vector normal a los nodos
    normnode = normal(geom);
else
    normnode = geom.normal;
end

if isfield(geom,'dsi') ~= 1
    % calcule el area superficial de los elementos
    [s,ds,dsi] = areas(geom);
else
    dsi = geom.dsi;
end

% calcule el centroide de la malla
if isfield(geom,'xc') ~= 1
    % calcule el centroide de la gota
   xc = centroide(geom); 
else
   xc = geom.xc; 
end

% referencie los nodos a la base ubicada en el centroide
nodesloc = geom.nodes - repmat(xc,[numnodes,1]);

integrand = inttensor(nodesloc,normnode);

% Calculo de la integral
inerttensor = (1/5).*inttrapeciomata(dsi,integrand);

return



%% Calculo del integrando

function tstnormal = inttensor(nodesxc,normnode)
numnodes = size(nodesxc,1);

xll = sum(nodesxc.*nodesxc,2);

% Calculo de los xgorro_i*xgorro_j
    dxx = nodesxc(:,1).^2;
    dxy = nodesxc(:,1).*nodesxc(:,2);
    dxz = nodesxc(:,1).*nodesxc(:,3);
    dyy = nodesxc(:,2).^2;
    dyz = nodesxc(:,2).*nodesxc(:,3);
    dzz = nodesxc(:,3).^2;

% Calculo de los terminos xixjxk

    tst2 = zeros(3,3,3,numnodes);
    
    tst2(1,1,1,:) = dxx.*nodesxc(:,1);
    tst2(1,1,2,:) = dxy.*nodesxc(:,1);
    tst2(2,1,1,:) = tst2(1,1,2,:);
    tst2(1,2,1,:) = tst2(1,1,2,:);
  
    tst2(1,1,3,:) = dxz.*nodesxc(:,1);
    tst2(3,1,1,:) = tst2(1,1,3,:);
    tst2(1,3,1,:) = tst2(1,1,3,:);
    
    tst2(2,1,2,:) = dyy.*nodesxc(:,1);
    tst2(2,2,1,:) = tst2(2,1,2,:);
    tst2(1,2,2,:) = tst2(2,1,2,:);
    
    tst2(2,1,3,:) = dyz.*nodesxc(:,1);
    tst2(3,1,2,:) = tst2(2,1,3,:);
    tst2(1,2,3,:) = tst2(2,1,3,:);
    tst2(1,3,2,:) = tst2(2,1,3,:);
    tst2(3,2,1,:) = tst2(2,1,3,:);
    tst2(2,3,1,:) = tst2(2,1,3,:);
   
    tst2(3,1,3,:) = dzz.*nodesxc(:,1);
    tst2(3,3,1,:) = tst2(3,1,3,:);
    tst2(1,3,3,:) = tst2(3,1,3,:);
    
    tst2(2,2,2,:) = dyy.*nodesxc(:,2);
    tst2(2,2,3,:) = dyz.*nodesxc(:,2);
    tst2(3,2,2,:) = tst2(2,2,3,:);
    
    tst2(3,2,3,:) = dzz.*nodesxc(:,2);
    tst2(3,3,2,:) = tst2(3,2,3,:);
    
    tst2(2,3,2,:) = dyy.*nodesxc(:,3); 
    tst2(2,3,3,:) = dyz.*nodesxc(:,3); 
    tst2(3,3,3,:) = dzz.*nodesxc(:,3);

% Calculo de los terminos Dijxk*(xlxl)

tst1 = zeros(3,3,3,numnodes);

    tst1(1,1,1,:) = nodesxc(:,1).*xll;
    tst1(1,1,2,:) = nodesxc(:,2).*xll;
    tst1(2,1,1,:) = 0;
    tst1(1,2,1,:) = 0;
  
    tst1(1,1,3,:) = nodesxc(:,3).*xll;
    tst1(3,1,1,:) = 0;
    tst1(1,3,1,:) = 0;
    
    tst1(2,1,2,:) = 0;
    tst1(2,2,1,:) = nodesxc(:,1).*xll;
    tst1(1,2,2,:) = 0;
    
    tst1(2,1,3,:) = 0;
    tst1(3,1,2,:) = 0;
    tst1(1,2,3,:) = 0;
    tst1(1,3,2,:) = 0;
    tst1(3,2,1,:) = 0;
    tst1(2,3,1,:) = 0;
   
    tst1(3,1,3,:) = 0;
    tst1(3,3,1,:) = nodesxc(:,1).*xll;
    tst1(1,3,3,:) = 0;
    
    tst1(2,2,2,:) = nodesxc(:,2).*xll;
    tst1(2,2,3,:) = nodesxc(:,3).*xll;
    tst1(3,2,2,:) = 0;
    
    tst1(3,2,3,:) = 0;
    tst1(3,3,2,:) = nodesxc(:,2).*xll;
    
    tst1(2,3,2,:) = 0; 
    tst1(2,3,3,:) = 0; 
    tst1(3,3,3,:) = nodesxc(:,3).*xll;

% Calculo de (tst1 - tst2)

tsttot = (tst1 - tst2);

% Ejecute tijk Nk para cada punto de integracion
    % Entrada tipo lista
    % Arreglo de valores de la normal para ejecutar producto en paralelo
    temp1 = permute(reshape(normnode',1,3*numnodes),[1 3 2]);
    normarray = reshape(repmat(temp1,[3 3 1]),[3 3 3 numnodes]);
    tstnormal = sum(tsttot.*normarray,3);

   % Permute para que quede array de 3 x 3 x NumQuadPoints*NumElements
    tstnormal = permute(tstnormal,[1 2 4 3]);
    
return
