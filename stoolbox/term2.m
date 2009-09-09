% calcula explicitamente el termino 2 de sulfactantes
function term2val = term2(gamma,u_s,struct)

if isfield(struct,'element2node') ~= 1
    % no hay tabla de elementos-nodos calculela
    ele2node = element2node(geom.elements);
else
    ele2node = struct.element2node;
end

if isfield(struct,'dsi') ~= 1
    % no hay tabla de elementos-nodos calculela
    [s,ds,dsi] = areas(struct);
else
    dsi = struct.dsi;
end

if isfield(struct,'normal') ~= 1 || isfield(struct,'normalele') ~= 1
    % calcule los vetores normals a los nodos y los elementos
    [normnode,normele] = normal(struct);
else
    normele = struct.normalele;
end

numnodes = size(struct.nodes,1);
gu_s = u_s.*repmat(gamma,[1 3]);

% Coordenadas nodo 4, 5 , 6 intermedios y 7 central
n4 = (struct.nodes(struct.elements(:,1),:) + struct.nodes(struct.elements(:,2),:)).*0.5;
n5 = (struct.nodes(struct.elements(:,3),:) + struct.nodes(struct.elements(:,2),:)).*0.5;
n6 = (struct.nodes(struct.elements(:,3),:) + struct.nodes(struct.elements(:,1),:)).*0.5;
n7 = struct.nodes(struct.elements(:,1),:).*(1/3) + ...
    struct.nodes(struct.elements(:,2),:).*(1/3) + ...
    struct.nodes(struct.elements(:,3),:).*(1/3);

% Diferencia entre nodo 7 e intermedios (genericamente n7_ni)
n7_n4 = n7 - n4;
n7_n5 = n7 - n5;
n7_n6 = n7 - n6;

n7_n4n = normesp(n7_n4);
n7_n5n = normesp(n7_n5);
n7_n6n = normesp(n7_n6);

% Productos cruz de los vectores respecto a sus normales (normalizados)
n74cnormal = crossv(n7_n4,normele);
n75cnormal = crossv(n7_n5,normele);
n76cnormal = crossv(n7_n6,normele);

n74cnormal = univect(n74cnormal);
n75cnormal = univect(n75cnormal);
n76cnormal = univect(n76cnormal);

% u_s en nodo 4,5 y 6 intermedios
gu_s4 = (gu_s(struct.elements(:,1),:) + gu_s(struct.elements(:,2),:)).*0.5;
gu_s5 = (gu_s(struct.elements(:,3),:) + gu_s(struct.elements(:,2),:)).*0.5;
gu_s6 = (gu_s(struct.elements(:,3),:) + gu_s(struct.elements(:,1),:)).*0.5;
gu_s7 = gu_s(struct.elements(:,1),:).*(1/3) + ...
    gu_s(struct.elements(:,2),:).*(1/3) + ...
    gu_s(struct.elements(:,3),:).*(1/3);

term2val = zeros(numnodes,1);
for k=1:numnodes
    % extraga los elementos vecinos al nodo la
    elek = ele2node{k};    
    inttot = 0;
    for z = 1:size(elek,2)
        % extraiga los nodos del elemento z
        nodz = struct.elements(elek(z),:);
        % encuentre la posicion del nodo k en el elemento elek(z)
        posz = find(nodz == k);      
        if posz == 1
        % nodo 1 local
            % vectores b1 y b2
            b1 = n74cnormal(elek(z),:);
            b2 = -n76cnormal(elek(z),:);
            m1 = n7_n4n(elek(z),:);
            m2 = n7_n6n(elek(z),:);
            intgus1 = 0.5.*m1.*(gu_s7(elek(z),:) + gu_s4(elek(z),:));
            intgus2 = 0.5.*m2.*(gu_s7(elek(z),:) + gu_s6(elek(z),:));
        elseif posz == 2
        % nodo 2 local
            % vectores b1 y b2
            b1 = n75cnormal(elek(z),:);
            b2 = -n74cnormal(elek(z),:);
            m1 = n7_n5n(elek(z),:);
            m2 = n7_n4n(elek(z),:);
            intgus1 = 0.5.*m1.*(gu_s7(elek(z),:) + gu_s5(elek(z),:));
            intgus2 = 0.5.*m2.*(gu_s7(elek(z),:) + gu_s4(elek(z),:));
        elseif posz == 3
        % nodo 3 local
            % vectores b1 y b2
            b1 = -n75cnormal(elek(z),:);
            b2 = n76cnormal(elek(z),:);
            m1 = n7_n5n(elek(z),:);
            m2 = n7_n6n(elek(z),:);
            intgus1 = 0.5.*m1.*(gu_s7(elek(z),:) + gu_s5(elek(z),:));
            intgus2 = 0.5.*m2.*(gu_s7(elek(z),:) + gu_s6(elek(z),:));
        end
        inttot = inttot + sum(b1.*intgus1,2) + sum(b2.*intgus2,2);
    end
    term2val(k) = inttot;    
end
term2val = term2val./dsi;
    
