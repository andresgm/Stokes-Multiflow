% Rutina de adapatacion de lowenberg and hinch
% TODO: refinar la rutina
function veladapt = meshadapt2(geom,adaptparms)

% recupere las propiedades geometricas
curvnodes = geom.curv;
normalnodes = geom.normal;
nodes = geom.nodes;
nodecon2node = geom.nodecon2node;
dsi = geom.dsi;
numnodes = geom.numnodes;
numelements = size(geom.elements,1);

% recupere las opciones y parametros de la adaptacion de malla
lamda = adaptparms.lamda;
% psi = adaptparms.psi;

curvop = abs(curvnodes).^1.5;

pj = projtensor(normalnodes);
veladapt = zeros(numnodes,3);
lmin = zeros(numnodes,1);
if geom.numdrops ~= 1
    % Adaptacion de lowenberg
    for j = 1:geom.numdrops
        nodestemp = geom.nodes;
        nodestemp(geom.nnodesdrop(j,1):geom.nnodesdrop(j,2),:) = 0;
        nodestemp(~any(nodestemp,2),:) = [];
        for k = geom.nnodesdrop(j,1):geom.nnodesdrop(j,2)
           % calcule para el punto respeto de los nodos de las otras gotas
           cantnodes = geom.numnodes - (geom.nnodesdrop(j,2)- geom.nnodesdrop(j,1) + 1);
           lmin(k) = min(normesp(repmat(geom.nodes(k,:),[cantnodes 1]) - ...
                     nodestemp));
%            if j == 1  
%                cantnodes = length(geom.nnodesdrop(2,1):geom.nnodesdrop(2,2));
%                lmin(k) = min(normesp(repmat(geom.nodes(k,:),[cantnodes 1]) - ...
%                geom.nodes(geom.nnodesdrop(2,1):geom.nnodesdrop(2,2),:)));  
%            elseif j == 2
%                cantnodes = length(geom.nnodesdrop(1,1):geom.nnodesdrop(1,2));
%                lmin(k) = min(normesp(repmat(geom.nodes(k,:),[cantnodes 1]) - ...
%                geom.nodes(geom.nnodesdrop(1,1):geom.nnodesdrop(1,2),:)));  
%            end
        end
    end

    for i = 1:numnodes
        % extraiga los nodos adyacentes al inode
        jnode = nodecon2node{i};
        deltaxji = nodes(jnode,:) - repmat(nodes(i,:),[size(jnode,1) 1]);
        ndeltaxji = normesp(deltaxji);
        hj = 1./lmin(jnode);
        dsij = dsi(jnode);
        % part1 = sum(1 + 4|xj-xi|/h(xj) + abs(k(xj))^1.5)*dsj*(xj - xi)
        part1 = sum(repmat((1 + 4.*ndeltaxji.*hj + curvop(jnode)).*dsij,[1 3]).*deltaxji,1);
        % multiplique por el tensor de proyeccion
        veladapt(i,:) = (numelements/geom.numdrops)^1.5/(300.*(1 + lamda)).*(pj(:,:,i)*part1')';
    end
else
    % Adaptacion para el caso cuando solo es una gota
    if strcmp(adaptparms.flow,'inf') == 1
        for i = 1:numnodes
            % extraiga los nodos adyacentes al inode
            jnode = nodecon2node{i};
            deltaxji = nodes(jnode,:) - repmat(nodes(i,:),[size(jnode,1) 1]);
            dsij = dsi(jnode);
            % part1 = sum(1 + 4|xj-xi|/h(xj) + abs(k(xj))^1.5)*dsj*(xj - xi)
            % adicionado constante 10 para curvatura
            part1 = sum(repmat((1 + 10.*curvop(jnode)).*dsij,[1 3]).*deltaxji,1);
            % multiplique por el tensor de proyeccion
            veladapt(i,:) = (numelements/geom.numdrops)^1.5/(300.*(1 + lamda)).*(pj(:,:,i)*part1')';
        end
    elseif strcmp(adaptparms.flow,'semiinf') == 1
        hjtot = 1./abs(nodes(:,3));
        for i = 1:numnodes
            % extraiga los nodos adyacentes al inode
            jnode = nodecon2node{i};
            deltaxji = nodes(jnode,:) - repmat(nodes(i,:),[size(jnode,1) 1]);
            ndeltaxji = normesp(deltaxji);
            hj = hjtot(jnode);
            dsij = dsi(jnode);
            % part1 = sum(1 + 4|xj-xi|/h(xj) + abs(k(xj))^1.5)*dsj*(xj - xi)
            part1 = sum(repmat((1 + 4.*ndeltaxji.*hj + 10.*curvop(jnode)).*dsij,[1 3]).*deltaxji,1);
            % multiplique por el tensor de proyeccion
            veladapt(i,:) = (numelements/geom.numdrops)^1.5/(300.*(1 + lamda)).*(pj(:,:,i)*part1')';
        end
        
%         for i = 1:numnodes
%             % extraiga los nodos adyacentes al inode
%             jnode = nodecon2node{i};
%             deltaxji = nodes(jnode,:) - repmat(nodes(i,:),[size(jnode,1) 1]);
%             dsij = dsi(jnode);
%             % part1 = sum(1 + 4|xj-xi|/h(xj) + abs(k(xj))^1.5)*dsj*(xj - xi)
%             part1 = sum(repmat((1 + curvop(jnode)).*dsij,[1 3]).*deltaxji,1);
%             % multiplique por el tensor de proyeccion
%             veladapt(i,:) = (numelements/geom.numdrops)^1.5/(300.*(1 + lamda)).*(pj(:,:,i)*part1')';
%         end
    end
end
