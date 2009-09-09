% Determina las coordenadas de los puntos de integracion de cada elemento
% 
% Entradas:
%      Elements: matriz de conectividad de los elementos dim(NE x 3) Linear
%           dim(NumElements x 6) Cuadratic
%      Nodes: Coordenadas de los nodos dim(NumNodes x 3)
%      NumQuadPoints: Numero de puntos de Integracion por elemento
%      IntType: Tipo de Integracion 1 Standard Gauss 2 Gauss Pozrikidis
% 
% Salidas:
%      QuadCoord: Coordenadas de los puntos de integracion en la base global
%           dim(NumElements x 3 x NumQuadPoints)
%      PhiZita: Matriz de funciones de forma evaluadas en cada punto de
%           integracion dim(3 x 9 x NumQuadPoints) Linear 
%           dim(3 x 18 x NumQuadPoints) Cuadratic
%      Weights: Vector de Pesos asociados a cada punto de integracion
%           dim(1 x NumQuadPoints)
%      PhiZitaWeights: matriz de Funciones de forma PhiZita multiplicada por 
%           cada peso para cada punto de integraci???n    
%           dim(3 x 9 x NumQuadPoints) Linear 
%           dim(3 x 18 x NumQuadPoints) Cuadratic
%      Zitas: Coordenadas locales Z1 y Z2 de cada punto de integracion
%           dim(2 x NumQuadPoints)

    
function [QuadCoord,PhiZita,Weights,PhiZitaWeights,Zitas] = QuadPoints(struct,numquadpoints,inttype)

NumQuadPoints = numquadpoints;

if nargin < 3
    % itnergracion estandar de gauss
IntType = 1;
else
    IntType = inttype;
end

Elements = struct.elements;
Nodes = struct.nodes;

% Numero de elementos
NumElements = size(Elements,1);
% Numero de nodos por elemento
NodesPerElement = size(Elements,2);

% Espacio en memoria
PhiZita = zeros(3,3*NodesPerElement,NumQuadPoints);
PhiZitaWeights = PhiZita;

if NodesPerElement == 3
    % Elemento Lineal
    GeomType = 1;
elseif NodesPerElement == 6
    % Elemento Cuadratico
    GeomType = 2;
end

% Coordenadas en la base isoparam???trica (Zita1,Zita2)* de los puntos de
% integraci???n y sus correspondientes pesos Weights
% * Zita3 = 1 - Zita1 - Zita2

if NumQuadPoints == 1
    % Tomado de Brebbia: Boundary Elements - An introductory course pg 289
    Zitas = [1/3;1/3];
    Weights = 1/2;
elseif NumQuadPoints == 3
    % Tomado de Brebbia: Boundary Elements - An introductory course pg 289
    Zitas = [0.5 0 0.5;0.5 0.5 0];
%     Zitas = [0 1 0;0 0 1];    
    Weights = [1/6 1/6 1/6];
elseif NumQuadPoints == 4
    % Tomado de Brebbia: Boundary Elements - An introductory course pg 289
    Zitas = [1/3 3/5 1/5 1/5;1/3 1/5 3/5 1/5]; 
    Weights = 0.5.*[-9/16 25/48 25/48 25/48];        
elseif NumQuadPoints == 7 && IntType == 1
    % 7 Puntos de integracion para Gauss Standard
    % Tomado de Brebbia: Boundary Elements - An introductory course pg 289
    % Valores con mas cifras decimales tomados de:
    % CRC The Finite Element Method Using Matlab pg 173
    Zitas = [1/3 0.7974269853531 0.1012865073235 0.1012865073235 0.0597158717898 0.4701420641051 0.4701420641051;...
        1/3 0.1012865073235 0.7974269853531 0.1012865073235 0.4701420641051 0.0597158717898 0.4701420641051];
    Weights = 0.5.*[0.225 0.1259391805448 0.1259391805448 0.1259391805448 ... 
        0.1323941527885 0.1323941527885 0.1323941527885];
elseif NumQuadPoints == 7 && IntType == 2
    % 7 Puntos de Integracion para Gauss Pozrikidis
    % Tomado de Pozrikidis Boundary Integral Formulation for Linearized
    % Viscous Flow pg 177
    Zitas = [1 0 0 0.5 0 0.5 1/3;0 1 0 0.5 0.5 0 1/3];
    Weights = 0.5.*[3/60 3/60 3/60 8/60 8/60 8/60 27/60];
elseif NumQuadPoints == 9
    % Tomado de C. Pozrikidis, Numerical Computation in Science and 
    % Engineering Chap Numerical Integration
    al = 0.124949503233232;
    qa = 0.165409927389841;
    rh = 0.797112651860071;
    de = 0.437525248383384;
    ru = 0.037477420750088;
    o1 = 0.205950504760887;
    o2 = 0.063691414286223;

    Zitas = [de de al ru qa qa rh rh ru;de al de qa ru rh qa ru rh];
    Weights = 0.5.*[o1 o1 o1 o2 o2 o2 o2 o2 o2];
end

% Calculo del valor de las funciones de forma en los puntos de integracion

% for i = 1:NumQuadPoints
%     PhiZita(:,:,i) = PhiShape(Zitas(:,i),GeomType);
%     PhiZitaWeights(:,:,i) = PhiZita(:,:,i).*Weights(i);
% end   

% Reorganice los nodos
nodo1 = Nodes(Elements(:,1),:);
nodo2 = Nodes(Elements(:,2),:);
nodo3 = Nodes(Elements(:,3),:);

QuadCoord = zeros(NumElements,3,NumQuadPoints);
for i = 1:NumQuadPoints
    % invoque interpolacion lineal
    
    QuadCoord(:,:,i) = lininterp(nodo1,nodo2,nodo3,Zitas(:,i));
end   

% otra opcion

% Calculo de las coordenadas de los puntos de integracion.
% Dim(i,j,k) = (NumElement x 3 x NuMQuadPoints)
% Cada slide i,j corresponde a las coordenadas del k-esimo punto de
% integracion de cada i-esimo elemento

% QuadCoord = PhiProduct(Nodes,PhiZita,Elements);


    
