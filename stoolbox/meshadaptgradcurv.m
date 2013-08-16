% Rutina de adaptaci?n de malla pasiva usando metodos de gradiente iterativos.
% Metodo basado en Zinchenco et al. 1997.
function veladapt = meshadaptgradcurv(geom,velnormal,veladapt,dt)

numnodes = geom.numnodes;

fi = zeros(numnodes,3);
maxit = 1000;
tolerance = 1e-6;
flag0 = 1;

Vel = velnormal + veladapt...
    - repmat(sum(veladapt.*geom.normal,2),[1 3]).*geom.normal;

fdev = FdeV(geom,Vel,dt);

for iter = 1:maxit
   fdevold = fdev;
   for i = 1:numnodes
      suma = zeros(1,3);
      for j = 1:size(geom.nodecon2node{i},1)
         xij = geom.nodes(geom.nodecon2node{i}(j),:) - geom.nodes(i,:);
         producto = (xij*(Vel(geom.nodecon2node{i}(j),:)-Vel(i,:))')*xij;
         suma = suma + producto;
      end
      fi(i,:) = (eye(3)-geom.normal(i,:)'*geom.normal(i,:))*suma';
   end

   eps = 1;
   flag = 1;
   [dfdeps,dfdepsprima] = dFdepsfast(geom,fi,Vel,eps);
   for k = 1:maxit
      eps = eps - dfdeps/dfdepsprima;
      [dfdeps,dfdepsprima] = dFdepsfast(geom,fi,Vel,eps);
      if (dfdeps < tolerance)
         flag = 0;
         break;
      end
   end

   if flag == 0
      Vel = Vel + eps*fi;
   else
      error('Newton-Raphson para calculo de Epsilon no convergio')
   end
   
   fdev = FdeV(geom,Vel,dt);
   if (fdevold - fdev < tolerance)
         flag0 = 0;
         break;
   end
      
end

if flag0 == 0
   veladapt = Vel - velnormal;
else
   error('Algoritmo de minimizacion no convergio')
end

end

function [dfdeps,dfdepsprima] = dFdepsfast(geom,f,Vel,eps)

   numvertices = size(geom.vertices,1);
   
   dfdeps = 0;
   dfdepsprima = 0;
   parfor i = 1:numvertices
      xij = geom.nodes(geom.vertices(i,2),:) - geom.nodes(geom.vertices(i,1),:);
      producto1 = xij*(Vel(geom.vertices(i,2),:)+eps*f(geom.vertices(i,2),:)...
         -Vel(geom.vertices(i,1),:)-eps*f(geom.vertices(i,1),:))';
      producto2 = xij*(f(geom.vertices(i,2),:)-f(geom.vertices(i,1),:))';
      dfdeps = dfdeps + producto1*producto2;
      dfdepsprima = dfdepsprima + producto2*producto2;
   end

end

% function dfdeps = dFdeps(geom,f,Vel,eps)
% 
%    numvertices = size(geom.vertices,1);
%    
%    suma = 0;
%    for i = 1:numvertices
%       xij = geom.nodes(geom.vertices(i,2),:) - geom.nodes(geom.vertices(i,1),:);
%       producto1 = xij*(Vel(geom.vertices(i,2),:)+eps*f(geom.vertices(i,2),:)...
%          -Vel(geom.vertices(i,1),:)-eps*f(geom.vertices(i,1),:))';
%       producto2 = xij*(f(geom.vertices(i,2),:)-f(geom.vertices(i,1),:))';
%       suma = suma + producto1*producto2;
%    end
%    dfdeps = suma;
% end      
% 
% function dfdepsprima = dFdepsprima(geom,f)
% 
%    numvertices = size(geom.vertices,1);
%    
%    suma = 0;
%    for i = 1:numvertices
%       xij = geom.nodes(geom.vertices(i,2),:) - geom.nodes(geom.vertices(i,1),:);
%       producto = xij*(f(geom.vertices(i,2),:)-f(geom.vertices(i,1),:))';
%       suma = suma + producto*producto;
%    end
%    dfdepsprima = suma;
% end

function fdev = FdeV(geom,Vel,dt)

    c1 = 0.25;
    c2 = 2;

   numvertices = size(geom.vertices,1);
   numnodes = geom.numnodes;
   
   k2max = zeros(numnodes,1);
   k2min = zeros(numnodes,1);
   phi = zeros(numnodes,1);
   
   for i = 1:numnodes
      k2 = [];
      for j = 1:size(geom.nodecon2node{i},1)
         xij = geom.nodes(geom.nodecon2node{i}(j),:) - geom.nodes(i,:);
         k2(j) = (geom.normal(geom.nodecon2node{i}(j),:)-geom.normal(i,:))*...
             (geom.normal(geom.nodecon2node{i}(j),:)-geom.normal(i,:))'...
             /(xij*xij');
      end
      k2max(i) = max(k2);
      k2min(i) = min(k2);
      phi(i) = k2max(i) - k2min(i)+ 1;
   end      
   
   geomvar = geom;
   geomvar.nodes = geomvar.nodes + dt*Vel;
   
   % calcule el vector normal a cada nodo
    normalandgeoopt.normal = 1;
    normalandgeoopt.areas = 1;
    normalandgeoopt.vol = 1;
    geomprop = normalandgeo(geom,normalandgeoopt,1);
    geomvar.normalele = geomprop.normalele;
    geomvar.normal = geomprop.normal;
    geomvar.dsi = geomprop.dsi;
    geomvar.ds = geomprop.ds;
    geomvar.s = geomprop.s;
    geomvar.vol = geomprop.vol;
    geomvar.jacmat = geomprop.jacmat;
    geomvar.g = geomprop.g;
    paropt.tipo = 'extended';
    [geomvar.curv,geomvar.normal,geomvar.Kg] = ...
        curvparaboloid(geomvar,paropt);
    dnormaldt = (geomvar.normal-geom.normal)./dt;
   
   term1 = 0;
   term2 = 0;
   for i = 1:numvertices
       % Primer termino 5.3 Cusping, capture, and breakup...
      xij = geom.nodes(geom.vertices(i,2),:) - ...
          geom.nodes(geom.vertices(i,1),:);
      der1 = (Vel(geom.vertices(i,2),:)-Vel(geom.vertices(i,1),:))*...
          (geom.normal(geom.vertices(i,2),:)-...
          geom.normal(geom.vertices(i,1),:))';
      der2 = xij*...
          (dnormaldt(geom.vertices(i,2),:)-dnormaldt(geom.vertices(i,1),:))';
      term1 = term1 + ((der1+der2)^2)/((xij*xij')^2);
      % Segundo termino 5.3 ...
      phii = phi(geom.vertices(i,1));
      phij = phi(geom.vertices(i,2));
      sumphis = (phii+phij)/((xij*xij')^2);
      dxijdt = ...
          4*(xij*(Vel(geom.vertices(i,2),:)-Vel(geom.vertices(i,1),:))')^2;
      term2 = term2 + sumphis*dxijdt;
   end
   
   for i = 1:geom.numelements
       % Tercer termino 5.3 ...
       phi1 = phi(geom.elements(i,1));
       phi2 = phi(geom.elements(i,2));
       phi3 = phi(geom.elements(i,3));
       phiave = (phi1 + phi2 + phi3)/3;
       cross
   
   fdev = term1 + c1*term2;
end