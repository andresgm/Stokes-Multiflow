% Calcula la matriz de laplace beltrami
% Calcula la curvatura Gaussiana mediante el oerador discreto de laplace
% beltrami
% TODO armar area voronio segun eq(3) Garimella-2003-curvature.pdf
function [l,Kg] = laplacebeltramimat(geom,type)

% Compute de laplace beltrami of a scalar field scfld for a triangulated elements
% Adapted from compute_geometric_laplacian function of Gabriel Peyre's
% graphtoolbox
% http://www.mathworks.com/matlabcentral/fileexchange ...
% ... /loadFile.do?objectId=5355&objectType=FIlE

% recupere las opciones

if nargin < 3
    %   TODO activar suo de area mixta y voronio
    type = 'voronio';
elseif strcmp(type,'voronio') == 1
    
end

if isfield(geom,'element2node') ~= 1
    % no se paso el cell de elementos a nodos calculelo
    ele2node = element2node(geom.elements);
else
    ele2node = geom.element2node;
end

%recupere las variables
nodes = geom.nodes;
elements = geom.elements;
numnodes = max(elements(:));
numelements = size(elements,1);

% conformal laplacian
l = sparse(numnodes,numnodes);
sumtheta = zeros(numnodes,1);

for i = 1:numnodes
    ele2nodev = ele2node{i};
    thetat = 0;
    for b = ele2nodev
        % b is elements adyacent to i node
        bf = elements(b,:);
        % define edges
        if bf(1) == i
            v = bf(2:3);
        elseif bf(2) == i
            v = bf([1 3]);
        elseif bf(3) == i
            v = bf(1:2);
        end
        j = v(1);
        k = v(2);
        vi = nodes(i,:);
        vj = nodes(j,:);
        vk = nodes(k,:);
        % angles
        alpha = angleedges(vk-vi,vk-vj);
        beta = angleedges(vj-vi,vj-vk);
        thetat = thetat + angleedges(vi-vk,vi-vj);
        % add weight
        l(i,j) = l(i,j) + cot(alpha);
        l(i,k) = l(i,k) + cot(beta);
    end
    % element angle at node i
    sumtheta(i) = thetat;
end

l = l - diag(sum(l,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the angles of each vertex of an element
angles2ele = zeros(numelements,3);
parfor ele=1:numelements
    % extraiga los nodos del j element
    bf = elements(ele,:);
       
    i = bf(1);
    j = bf(2);
    k = bf(3);
    vi = nodes(i,:);
    vj = nodes(j,:);
    vk = nodes(k,:);

    % angles
    a1 = angleedges(vi-vj,vi-vk);
    a2 = angleedges(vj-vi,vj-vk);
    a3 = angleedges(vk-vi,vk-vj);

    angles2ele(ele,:) = [a1 a2 a3];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute circuncentric areas for each node
area2nodo = area2nodes(nodes,elements,angles2ele);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Baricentric area only
% ZitaVectnodes = [1/3;1/3];
% Weights = 1/2;
% [G,g1,g2,g3,dX1Z1,dX1Z2,dX2Z1,dX2Z2,dX3Z1,dX3Z2] = MetricTrans(nodes,elements,ZitaVectnodes);
% % Area superficial de cada Elemento y Total
% dS = G.*Weights;
% % Asigne la 3ra parte a cada nodo y sume
% dSi = zeros(numnodes,1);
% for i=1:numnodes
%    % Extraiga los elementos del i nodo
%    iEle = ele2node{i};
%    dSi(i) = (1/3)*sum(dS(iEle),1);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parfor j = 1:numnodes
   l(:,j) = l(:,j)./(2.*area2nodo);
end

% Gaussian Curvature
Kg = ((2.*pi) - sumtheta)./area2nodo;

%% Auxiliar functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute angle between 2 edges
function teta = angleedges(u,v)

du = sqrt(sum(u.^2));
dv = sqrt(sum(v.^2));
du = max(du,eps); dv = max(dv,eps);
% cos(teta) = <u,v>/(|u|*|b|)
teta = acos(sum(u.*v)/(du*dv));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute de area surrounding a node for a circuncentric dual
function area2nodo = area2nodes(nodes,elements,angles2ele)

numelements = size(elements,1);
numnodes = size(nodes,1);

% Coordenadas nodos 1, 2 y 3 vertices cada elemento
nodo1 = nodes(elements(:,1),:);
nodo2 = nodes(elements(:,2),:);
nodo3 = nodes(elements(:,3),:);
% Coordenadas nodo 4, 5 y 6 intermedios
nodo4 = (nodes(elements(:,1),:) + nodes(elements(:,2),:)).*0.5;
nodo5 = (nodes(elements(:,3),:) + nodes(elements(:,2),:)).*0.5;
nodo6 = (nodes(elements(:,3),:) + nodes(elements(:,1),:)).*0.5;
% circuncentro de cada triangulo 
nodo7 = circuncentre(nodo1,nodo2,nodo3);

% Subareas a cada nodo
sarea1 = trianglearea(nodo1,nodo4,nodo7);
sarea2 = trianglearea(nodo4,nodo2,nodo7);
sarea3 = trianglearea(nodo2,nodo5,nodo7);
sarea4 = trianglearea(nodo7,nodo5,nodo3);
sarea5 = trianglearea(nodo7,nodo3,nodo6);
sarea6 = trianglearea(nodo7,nodo6,nodo1);

area1 = sarea1 + sarea6;
area2 = sarea2 + sarea3;
area3 = sarea4 + sarea5;

% Indexing areas for each local node

area2nodolocal = [area1 area2 area3];

% Area of each element
elementarea = trianglearea(nodo1,nodo2,nodo3);

area2nodo = zeros(numnodes,1);

for i=1:numelements
    nodesele = elements(i,:);
    anglesele = angles2ele(i,:);
    % pregunte si algun angulo del elemento es obtuso
    obtind = anglesele>pi/2;
    if anglesele(obtind) > pi/2
        % el triangulo es obtuso
        % para el nodo del triangulo obstuso        
        area2nodo(nodesele(obtind)) = area2nodo(nodesele(obtind)) + elementarea(i)./2;
        % para los nodos del triangulo no obtuso
        area2nodo(nodesele(obtind == 0)) = area2nodo(nodesele(obtind == 0)) + elementarea(i)./4;
    else
        % el triangulo no es obtuso Voronoi area safe
        area2nodo(nodesele) = area2nodo(nodesele) + area2nodolocal(i,:)';
    end

end
