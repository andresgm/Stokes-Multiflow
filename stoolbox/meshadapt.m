function veladapt = meshadapt(geom,adaptparms,velatnode)

% recupere las propiedades geometricas
curvatnodes = geom.curv;
normalatnodes = geom.normal;
nodes = geom.nodes;
elements = geom.elements;
elesingularIndex = geom.element2node;
nodecon2node = geom.nodecon2node;
dsi = geom.dsi;
ds = geom.ds;
numelements = size(geom.elements,1);
edgeIndex = geom.edgeindex;
numnodes = geom.numnodes;

% recupere las opciones y parametros de la adaptacion de malla
Lmin = adaptparms.lmin;
Lamda = adaptparms.lamda;
psi = adaptparms.psi;
adpopt = adaptparms.opt; 

if adpopt == 0
    % no adaptacion de malla
    veladapt = zeros(numnodes,3);
elseif adpopt == 1
    % adaptacion de malla segun Zinchenko 1997
%     psi = 0.2;
    prtenatnode = projtensor(normalatnodes);
    veladapt = zeros(numnodes,3);
    for k = 1:numnodes
        % extraiga los nodos vecinos al inode
        nodeadj = nodecon2node{k};
        cte1 = (1 + 4.*normesp(nodes(nodeadj,:) - repmat(nodes(k,:),[size(nodeadj,1) 1]))./Lmin...
            + abs((2.*curvatnodes(nodeadj)).^1.5)).*(dsi(nodeadj));
        vect = nodes(nodeadj,:) - repmat(nodes(k,:),[size(nodeadj,1) 1]);
        velad = (psi*numelements^1.5/(300*(1 + Lamda))).*prtenatnode(:,:,k)*(sum(repmat(cte1,[1 3]).*vect,1))';
        veladapt(k,:) = velad';
    end
elseif adpopt == 2
    Factor = 2;
% adaptacion de malla basada en el centroide
    mc = zeros(numnodes,3);
    for i=1:numnodes
       % extraiga los elementos del i nodo
       iele = elesingularIndex{i};
       temp1 = zeros(1,3);
       for k=1:size(iele,1)
          nodeele = elements(iele(k),:);
          temp1 = sum(nodes(nodeele,:),1)./3.*ds(iele(k)) + temp1;
       end
       % centro de masa del i nodo
       mc(i,:) = temp1/sum(ds(iele,1));
    end
    mc_nodes = mc - nodes;
    prtenatnode = projtensor(normalatnodes);
    temp1 = sum(prtenatnode.*repmat(permute(permute(mc_nodes',[1 3 2]),[2 1 3]),[3 1 1]),2);
    mc_nodest = reshape(temp1,3,numnodes)';
    veladapt = Factor.*mc_nodest;
elseif adpopt == 3
% Determine las longitudes de los bordes
    numedges = size(edgeIndex,1);
    l = zeros(numedges,1);
    for k=1:numedges
        % extraiga los nodos del kedge
        kedge = edgeIndex(k,:);
        % Longitud del kedge
        l(k) = normesp(nodes(kedge(1),:) - nodes(kedge(2),:));
    end
    Ktec = 1;
    Kte = 1;
    for i=1:numnodes
        edgeInd = find(edgeIndex(:,1) == i);
        h(i) = mean(l(edgeIndex(edgeInd)));        
        B(i) = Kte*h(i);
    end
    
    prtenatnode = projtensor(normalatnodes);
    veladapt = zeros(numnodes,3);
    for k = 1:numnodes
        % extraiga los nodos vecinos al inode
        nodeadj = nodecon2node{k};
        cte1 = (1 + B(nodeadj)' + Ktec.*abs(curvatnodes(nodeadj)).^1.5).*(dsi(nodeadj));
        vect = nodes(nodeadj,:) - repmat(nodes(k,:),[size(nodeadj,1) 1]);
        avvel = sum(velatnode(nodeadj,:),1)./size(nodeadj,1);
        velad = prtenatnode(:,:,k)*(sum(repmat(cte1,[1 3]).*vect,1) - velatnode(k,:) + avvel)';
        veladapt(k,:) = velad';
    end
elseif adpopt == 4
    % Bazhlekov non singular contour integral
     % adaptacion de malla segun Zinchenko 1997
%     psi = 0.2;
    prtenatnode = projtensor(normalatnodes);
    veladapt = zeros(numnodes,3);
    a = 1;
    b = 4;
    c = 10;
    for k = 1:numnodes
        % extraiga los nodos vecinos al inode
        nodeadj = nodecon2node{k};
        % calcule las distancias minimas del inode al jnode
        hj = normesp(nodes(nodeadj,:) - repmat(nodes(k,:),[size(nodeadj,1) 1]));
        cte1 = (a + b.*hj + c.*abs((2.*curvatnodes(nodeadj)).^1.5)).*(dsi(nodeadj));
        vect = nodes(nodeadj,:) - repmat(nodes(k,:),[size(nodeadj,1) 1]);
        % velocidad promedio us
        Us = sum([velatnode(k,:);velatnode(nodeadj,:)],1)./(size(nodeadj,1) + 1);
        temp1 = (sum(repmat(cte1,[1 3]).*vect,1) - velatnode(k,:) + Us);
        velad = prtenatnode(:,:,k)*(temp1');
        veladapt(k,:) = velad';
    end
elseif adpopt == 5
    % Bazhlekov non singular contour integral
     % adaptacion de malla segun Zinchenko 1997
%     psi = 0.2;
    prtenatnode = projtensor(normalatnodes);
    veladapt = zeros(numnodes,3);
    a = 1;
    b = 4;
    c = 1;
    parfor k = 1:numnodes
        % extraiga los nodos vecinos al inode
        nodeadj = nodecon2node{k};
        % h de cada nodo a la pared
        hj = nodes(nodeadj,3);
%         cte1 = (a + b.*normesp(nodes(nodeadj,:) - repmat(nodes(k,:),[size(nodeadj,1) 1]))./hj...
%             + c.*abs((2.*curvatnodes(nodeadj)).^1.5)).*(dsi(nodeadj));
        cte1 = (a + b./hj + c.*abs((2.*curvatnodes(nodeadj)).^1.5)).*(dsi(nodeadj));
        vect = nodes(nodeadj,:) - repmat(nodes(k,:),[size(nodeadj,1) 1]);
        % velocidad promedio us
        Us = sum([velatnode(k,:);velatnode(nodeadj,:)],1)./(size(nodeadj,1) + 1);
        temp1 = (sum(repmat(cte1,[1 3]).*vect,1) - velatnode(k,:) + Us);
        velad = prtenatnode(:,:,k)*(temp1');
        veladapt(k,:) = velad';
    end    
end
