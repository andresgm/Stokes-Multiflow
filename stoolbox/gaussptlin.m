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

    
function [quadcoord,phizita,weights,phizitaweights,zitas] = gaussptlin(struct,numquadpoints,inttype)


if nargin < 3
    % itnergracion estandar de gauss
    inttype = 1;
end

elements = struct.elements;
nodes = struct.nodes;

% numero de elementos
numelements = size(elements,1);
% numero de nodos por elemento
nodesperelement = size(elements,2);

% espacio en memoria
phizita = zeros(3,3*nodesperelement,numquadpoints);
phizitaweights = phizita;

% coordenadas en la base isoparam???trica (zita1,zita2)* de los puntos de
% integraci???n y sus correspondientes pesos weights
% * zita3 = 1 - zita1 - zita2

if numquadpoints == 1
    % tomado de Brebbia: Boundary elements - An introductory course pg 289
    zitas = [1/3;1/3];
    weights = 1/2;
elseif numquadpoints == 3
    % tomado de Brebbia: Boundary elements - An introductory course pg 289
    zitas = [0.5 0 0.5;0.5 0.5 0];
%     zitas = [0 1 0;0 0 1];    
    weights = [1/6 1/6 1/6];
elseif numquadpoints == 4
    % tomado de Brebbia: Boundary elements - An introductory course pg 289
    zitas = [1/3 3/5 1/5 1/5;1/3 1/5 3/5 1/5]; 
    weights = 0.5.*[-9/16 25/48 25/48 25/48];        
elseif numquadpoints == 7 && inttype == 1
    % 7 puntos de integracion para gauss Standard
    % tomado de Brebbia: Boundary elements - An introductory course pg 289
    % Valores con mas cifras decimales tomados de:
    % cRc the Finite element Method Using Matlab pg 173
    zitas = [1/3 0.7974269853531 0.1012865073235 0.1012865073235 0.0597158717898 0.4701420641051 0.4701420641051;...
        1/3 0.1012865073235 0.7974269853531 0.1012865073235 0.4701420641051 0.0597158717898 0.4701420641051];
    weights = 0.5.*[0.225 0.1259391805448 0.1259391805448 0.1259391805448 ... 
        0.1323941527885 0.1323941527885 0.1323941527885];
elseif numquadpoints == 7 && inttype == 2
    % 7 puntos de integracion para gauss pozrikidis
    % tomado de pozrikidis Boundary integral Formulation for Linearized
    % Viscous Flow pg 177
    zitas = [1 0 0 0.5 0 0.5 1/3;0 1 0 0.5 0.5 0 1/3];
    weights = 0.5.*[3/60 3/60 3/60 8/60 8/60 8/60 27/60];
elseif numquadpoints == 9
    % tomado de c. pozrikidis, numerical computation in Science and 
    % engineering chap numerical integration
    al = 0.124949503233232;
    qa = 0.165409927389841;
    rh = 0.797112651860071;
    de = 0.437525248383384;
    ru = 0.037477420750088;
    o1 = 0.205950504760887;
    o2 = 0.063691414286223;

    zitas = [de de al ru qa qa rh rh ru;de al de qa ru rh qa ru rh];
    weights = 0.5.*[o1 o1 o1 o2 o2 o2 o2 o2 o2];
end

% Reorganice los nodos
nodo1 = nodes(elements(:,1),:);
nodo2 = nodes(elements(:,2),:);
nodo3 = nodes(elements(:,3),:);

quadcoord = zeros(numelements,3,numquadpoints);
for i = 1:numquadpoints
    % invoque interpolacion lineal    
    quadcoord(:,:,i) = lininterp(nodo1,nodo2,nodo3,zitas(:,i));
end   

