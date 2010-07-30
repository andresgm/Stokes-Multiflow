% Calcula la metrica de transformacion para un conjunto de puntos del
% elemento ZitaVect en la base isoparametrica. Lo anterior para cada
% elemento de un arreglo de elementos Elements
% Tomado de Brebbia: Boundary elements - An introductory course pg 120-122

function [jacomp,dx1,dx2,dx3] = metrictrans(geom,zitavect,opcion)

% recupere las variables
numelements = size(geom.elements,1);
nodesperelement = size(geom.elements,2);
numpoints = size(zitavect,2);

elements = geom.elements;
nodes = geom.nodes;

% Derivadas locales de cada punto de integracion.
    dx1z1 = zeros(numelements,numpoints);
    dx1z2 = zeros(numelements,numpoints);
    dx2z1 = zeros(numelements,numpoints);
    dx2z2 = zeros(numelements,numpoints);
    dx3z1 = zeros(numelements,numpoints);
    dx3z2 = zeros(numelements,numpoints);
    
if nodesperelement == 3
    % caso de elementos lineales triangulares 3 nodos
    
    % Rutina calculo de las derivadas de x en el elemento isoparametrico
    elementstemp = reshape(elements',numelements*nodesperelement,1);
    x1 = nodes(elementstemp,1);
    x2 = nodes(elementstemp,2);
    x3 = nodes(elementstemp,3);
    x1 = reshape(x1,nodesperelement,numelements);
    x2 = reshape(x2,nodesperelement,numelements);
    x3 = reshape(x3,nodesperelement,numelements);
    
    % calculo de las derivadas dx/zi 
    [dx1z1, dx1z2, dx2z1, dx2z2, dx3z1, dx3z2,dx1z3,dx2z3,dx3z3] = ...
        dfzitalin(x1,x2,x3);
       
    % Replique para los puntos de integracion
    dx1z1 = repmat(dx1z1',1,numpoints);
    dx1z2 = repmat(dx1z2',1,numpoints);
    dx2z1 = repmat(dx2z1',1,numpoints);
    dx2z2 = repmat(dx2z2',1,numpoints);
    dx3z1 = repmat(dx3z1',1,numpoints);
    dx3z2 = repmat(dx3z2',1,numpoints);
    
    dx1z3 = repmat(dx1z3',1,numpoints);
    dx2z3 = repmat(dx2z3',1,numpoints);
    dx3z3 = repmat(dx3z3',1,numpoints);
    
    dx1.dz1 = dx1z1;
    dx1.dz2 = dx1z2;
    dx2.dz1 = dx2z1;
    dx2.dz2 = dx2z2;
    dx3.dz1 = dx3z1;
    dx3.dz2 = dx3z2;
    
    dx1.dz3 = dx1z3;
    dx2.dz3 = dx2z3;
    dx3.dz3 = dx3z3;
elseif nodesperelement == 6
    % elementos triangulares de segundo orden 6 nodos
 
    % nUevA RUtInA tenIenDO en cUentA QUe nODOS InteRmeDIOS nO cOIncIDen
    % cOn zItA 0.5
    [dx1z1, dx1z2, dx2z1, dx2z2, dx3z1, dx3z2] = ...
        dfzitacurv2(nodes,elements,zitavect,nodes);
end

% componentes (g1,g2,g3) de la normal
% Brebbia: Boundary elements - An introductory course pg 122
jacomp.g1 = dx2z1.*dx3z2 - dx2z2.*dx3z1;
jacomp.g2 = dx3z1.*dx1z2 - dx1z1.*dx3z2;
jacomp.g3 = dx1z1.*dx2z2 - dx2z1.*dx1z2;

% Jacobiano reducido o metrica de transformacion
% Brebbia: Boundary elements - An introductory course pg 122
jacomp.g = (jacomp.g1.^2 + jacomp.g2.^2 + jacomp.g3.^2).^0.5;

if nargin == 3
    % Jacobiano completo
    jacinv = zeros(3,3,numelements);
    parfor i = 1:numelements
       % extraiga y ensamble las componentes del Jacobiano de cada elemento
       jacobiano = [dx1.dz1(i) dx2.dz1(i) dx3.dz1(i);...
                dx1.dz2(i) dx2.dz2(i) dx3.dz2(i);...
                dx1.dz3(i) dx2.dz3(i) dx3.dz3(i)];
       % genere su inversa
        jacinv(:,:,i) = inv(jacobiano);
    end
    jacomp.jacinv = jacinv;

else
    jacomp.jacinv = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     % RUtInA ORIgInAl elementOS tRIAngUlAReS De SegUnDO ORDen 6 nODOS
%     for i = 1:numelements
%         % nodos que pertenecen al i-esimo elemento
%             nodeselement = elements(i,:);
%         % coordenadas x1, x2, x3 (Dim,6x1) de los nodos (1er termino>1er nodo)
%             x1 = nodes(nodeselement,1);
%             x2 = nodes(nodeselement,2);
%             x3 = nodes(nodeselement,3);
%         % Derivadas de las coordenadas globales (x1,x2,x3) respecto a las
%         % locales (z1,z2) (elementos triangulares curvilineos)
%         % calculo de las derivadas dx/zi
%             [dx1z1(i,:), dx1z2(i,:), dx2z1(i,:), dx2z2(i,:), dx3z1(i,:), dx3z2(i,:)] = dfzitacurv(x1,x2,x3,zitavect);
%     end
