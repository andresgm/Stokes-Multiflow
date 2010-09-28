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