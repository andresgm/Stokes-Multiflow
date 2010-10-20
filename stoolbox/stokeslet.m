%% Calcula la funcion de green en espacio libre para un punto x0 debido a
%% todos los punto x.
%%          entradas:
%%              x0: Vector de 1 x 3 de coordenadas del polo
%%              x: Matriz de n x 3 con n el numero de puntos. Cada fila es 
%%              la coordenada del i-esimo punto.
%%          salidas:
%%              g: es una hipermatriz de 3 x 3 x n donde para cada n se 
%%              obtiene la matriz de funciones de green de x0 debido a cada
%%              i-esimo punto.
%% g(x,x0) = 1/norm(x^) - xi^xj^/(norm(x^)^3) con
%% x^ = x - x0

%% Tomado de Pozrikidis Boundary Integral and singularity Methods for
%% Linearized Viscous Flow Chap 2 pg 22
%% Adaptado de codigo stokes Flow Pozrikidis - BeMLIB

function [g,closenode] = stokeslet(x0,x)

% numero de Puntos x
    numnodes = size(x,1);

% Asignacion de espacio hipermatriz de funciones de green
    g = zeros(3,3,numnodes);

% Vector x-x0
    xgorro = x - repmat(x0,numnodes,1);

% norma de x-x0 para cada componente y norma al cubo
    r = normesp(xgorro);
    % determina el indice del punto mas cercano a x0
    closenode = find(r == min(r) & r > eps);
    if size(closenode,1) ~= 1 && size(closenode,1) ~= 0
        closenode = closenode(1);
    elseif  size(closenode,1) == 0
        % el nodo mas cercano esta en el mismo lugar que x0
        closenode = find(r == min(r));
    end
    
    rcub = r.^3;

% Quitar Warning de division por cero
    warning off MATLAB:divideByZero

% Calculo de los xgorro_i*xgorro_j
    dxx = xgorro(:,1).^2;
    dxy = xgorro(:,1).*xgorro(:,2);
    dxz = xgorro(:,1).*xgorro(:,3);
    dyy = xgorro(:,2).^2;
    dyz = xgorro(:,2).*xgorro(:,3);
    dzz = xgorro(:,3).^2;

% Inverso de la distancia al cubo
    ri3 = 1.0./rcub(:);

% Calculo de los elementos de la diagonal del stokeslet
    g(1,1,:) = 1./r(:) + dxx.*ri3;
    g(2,2,:) = 1./r(:) + dyy.*ri3;
    g(3,3,:) = 1./r(:) + dzz.*ri3;
      
% Calculo de los elementos de la diagonal superior
    g(1,2,:) = dxy.*ri3;
    g(1,3,:) = dxz.*ri3;
    g(2,3,:) = dyz.*ri3;
      
% elementos Diagonal Inferior = Diagonal superior
    g(2,1,:) = g(1,2,:);
    g(3,1,:) = g(1,3,:);
    g(3,2,:) = g(2,3,:);
    
% Multiplique por -1/8pi

g = (-1/(8*pi)).*g;
