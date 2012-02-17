%% Calcula el tensor de proyeccion:
%% Pg = [I] - nn %%% Pg_ij = dij - ni*nj
%% Ecuacion 5.6.4 Pozrikidis: Boundary Integral and Singularity Methods for
%% Linearized Viscous Flow pg 157. Tambi?n consultar en
%% Cristini V. An Adaptive Mesh Algorithm for Evolving Surfaces Ecuaci?n 9
%% pg 447.
%%      Entradas:
%%          Normal: Array de Normales Dim(NumNodes,3)
%%      Salidas:
%%          Pg: Tensor de Proyeccion Dim(3,3,NumNodes);

function pg = projtensor(normal)

numnodes = size(normal,1);

% pg = zeros (3,3,numnodes);

normal = reshape(normal',1,3,numnodes);

nn(1,1,:) = normal(1,1,:).*normal(1,1,:);
nn(1,2,:) = normal(1,1,:).*normal(1,2,:);
nn(1,3,:) = normal(1,1,:).*normal(1,3,:);
nn(2,1,:) = normal(1,2,:).*normal(1,1,:);
nn(2,2,:) = normal(1,2,:).*normal(1,2,:);
nn(2,3,:) = normal(1,2,:).*normal(1,3,:);
nn(3,1,:) = normal(1,3,:).*normal(1,1,:);
nn(3,2,:) = normal(1,3,:).*normal(1,2,:);
nn(3,3,:) = normal(1,3,:).*normal(1,3,:);

i = eye(3);
i = repmat(i,[1 1 numnodes]);

pg = i - nn;
