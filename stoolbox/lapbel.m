% Calcula el operador discreto de Laplace-Beltrami de la funcion f calculado
% sobre la superficie discreta de elementos triangulares descritos por
% geom.nodes y geom.elements
function lapbemf = lapbel(geom,f)

% Entradas: geom contiene por lo menos geom.nodes con las coordenadas de los
% nodos y la matriz de coordinacion geom.elements, f funcion escalar definida
% sobre la superficie.
% Salida: valor del Laplace-Beltami de f sobre la superficie.

nodecon2node = geom.nodecon2node;
nodes = geom.nodes;
elements = geom.elements;
numnodes = max(elements(:));
lapbemf = zeros(numnodes,1);

for i = 1:numnodes
   jnb = nodecon2node{i};
   valencia = size(jnb,1);
   alpha = zeros(valencia,1);
   beta = zeros(valencia,1);
   avori = 0;
   lapbemi = 0;
   for j = 1:valencia
      if j == 1
         alpha(1) = ...
            angleedges(nodes(i,:)-nodes(jnb(end),:),nodes(jnb(j),:)-nodes(jnb(end),:));
      else
         alpha(j) = ...
            angleedges(nodes(i,:)-nodes(jnb(j-1),:),nodes(jnb(j),:)-nodes(jnb(j-1),:));
      end
      if j == valencia
         beta(valencia) = ...
            angleedges(nodes(i,:)-nodes(jnb(1),:),nodes(jnb(j),:)-nodes(jnb(1),:));
      else
         beta(j) = ...
            angleedges(nodes(i,:)-nodes(jnb(j+1),:),nodes(jnb(j),:)-nodes(jnb(j+1),:));
      end
      normxij = normesp(nodes(i,:)-nodes(jnb(j),:));
      avori = avori + (cot(alpha(j))+cot(beta(j)))*normxij^2;
      lapbemi = lapbemi + (cot(alpha(j))+cot(beta(j)))*(f(i)-f(jnb(j)));
   end
   avori = avori/8;
   lapbemf(i) = lapbemi/(2*avori);
end

end

function teta = angleedges(u,v)

   du = sqrt(sum(u.^2));
   dv = sqrt(sum(v.^2));
   du = max(du,eps); dv = max(dv,eps);
   % cos(teta) = <u,v>/(|u|*|b|)
   teta = acos(sum(u.*v)/(du*dv));
   
end