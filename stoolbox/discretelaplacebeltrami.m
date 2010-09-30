% Calcula la matriz del operador discreto de Laplace-Beltrami

function [L,Kg] = discretelaplacebeltrami(geom)

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
    dist2 = zeros(numnodes,1);
    for b = ele2nodev
        % b is elements adyacent to i node
        bf = elements(b,:);
        % define edges
        if bf(1) == i
            v = bf(2:3);
        elseif bf(2) == i
            v = bf([3 1]);
        elseif bf(3) == i
            v = bf(1:2);
        end
        j = v(1);
        k = v(2);
        vi = nodes(i,:);
        vj = nodes(j,:);
        vk = nodes(k,:);
        % angles
        beta = angleedges(vi-vk,vj-vk);
        alpha = angleedges(vi-vj,vk-vj);
        thetat = thetat + angleedges(vk-vi,vj-vi);
        dist2(j) = normesp(vi-vj)^2;
        % add weight
        l(i,j) = l(i,j) + cot(beta);
        l(i,k) = l(i,k) + cot(alpha);
    end
    % element angle at node i
    sumtheta(i) = thetat;
    
    Avor(i) = l(i,:)*dist2/8;
    
end

L = l - diag(sum(l,2));

for i = 1:numnodes
   L(i,:) = L(i,:)./(2.*Avor(i));
end

% Gaussian Curvature
Kg = ((2.*pi) - sumtheta)./Avor';



function teta = angleedges(u,v)

du = sqrt(sum(u.^2));
dv = sqrt(sum(v.^2));
du = max(du,eps); dv = max(dv,eps);
% cos(teta) = <u,v>/(|u|*|b|)
teta = acos(sum(u.*v)/(du*dv));