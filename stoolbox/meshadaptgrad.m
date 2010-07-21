% Rutina de adaptaci?n de malla pasiva usando metodos de gradiente iterativos.
% Metodo basado en Zinchenco et al. 1997.
function veladapt = meshadaptgrad(geom,velnormal,veladapt)

numnodes = geom.numnodes;

fi = zeros(numnodes,3);
maxit = 1000;
tolerance = 1e-6;
flag0 = 1;

Vel = velnormal + veladapt;

fdev = FdeV(geom,Vel);

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
   dfdeps = dFdeps(geom,fi,Vel,eps);
   dfdepsprima = dFdepsprima(geom,fi);
   for k = 1:maxit
      eps = eps - dfdeps/dfdepsprima;
      dfdeps = dFdeps(geom,fi,Vel,eps);
      dfdepsprima = dFdepsprima(geom,fi);
      if (dfdeps < tolerance)
         flag = 0;
         break;
      end
   end

   if flag == 0
      Vel = Vel + eps*fi;
   else
      error('Newton-Rapson para calculo de Epsilon no convergio')
   end
   
   fdev = FdeV(geom,Vel);
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

function dfdeps = dFdeps(geom,f,Vel,eps)

   numvertices = size(geom.vertices,1);
   
   suma = 0;
   for i = 1:numvertices
      xij = geom.nodes(geom.vertices(i,2),:) - geom.nodes(geom.vertices(i,1),:);
      producto1 = xij*(Vel(geom.vertices(i,2),:)+eps*f(geom.vertices(i,2),:)...
         -Vel(geom.vertices(i,1),:)-eps*f(geom.vertices(i,1),:))';
      producto2 = xij*(f(geom.vertices(i,2),:)-f(geom.vertices(i,1),:))';
      suma = suma + producto1*producto2;
   end
   dfdeps = suma;
end      

function dfdepsprima = dFdepsprima(geom,f)

   numvertices = size(geom.vertices,1);
   
   suma = 0;
   for i = 1:numvertices
      xij = geom.nodes(geom.vertices(i,2),:) - geom.nodes(geom.vertices(i,1),:);
      producto = xij*(f(geom.vertices(i,2),:)-f(geom.vertices(i,1),:))';
      suma = suma + producto*producto;
   end
   dfdepsprima = suma;
end

function fdev = FdeV(geom,Vel)

   numvertices = size(geom.vertices,1);
   
   suma = 0;
   for i = 1:numvertices
      xij = geom.nodes(geom.vertices(i,2),:) - geom.nodes(geom.vertices(i,1),:);
      producto = xij*(Vel(geom.vertices(i,2),:)-Vel(geom.vertices(i,1),:))';
      suma = suma + producto*producto;
   end
   fdev = 4*suma;
end