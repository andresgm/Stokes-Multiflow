% calcula la matriz del termino 2 de sulfactantes
function matterm2 = c_term2mat(struct,u_s)

if isfield(struct,'normal') ~= 1 || isfield(struct,'normalele') ~= 1 || ...
        isfield(struct,'dsi') ~= 1 ||  isfield(struct,'element2node') ~= 1
    % no esta calculado el vector normal a los elementos, calculelo
    jacomp = metrictrans(struct,[1/3;1/3]);
    [normnode,normele] = normal(struct,jacomp);
    [s,ds,dsi] = areas(struct);
    ele2node = element2node(struct.elements);
else
    normele = struct.normalele;
    dsi = struct.dsi;
    ele2node = struct.element2node;
end

numnodes = size(struct.nodes,1);

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
u_s4 = (u_s(struct.elements(:,1),:) + u_s(struct.elements(:,2),:)).*0.5;
u_s5 = (u_s(struct.elements(:,3),:) + u_s(struct.elements(:,2),:)).*0.5;
u_s6 = (u_s(struct.elements(:,3),:) + u_s(struct.elements(:,1),:)).*0.5;
u_s7 = u_s(struct.elements(:,1),:).*(1/3) + ...
    u_s(struct.elements(:,2),:).*(1/3) + ...
    u_s(struct.elements(:,3),:).*(1/3);

%%% Cambia a la vectorizada
matterm2 = zeros(numnodes,numnodes);
for k=1:numnodes
    % extraiga los elementos del nodo k
    elek = ele2node{k};
       
    for z =1:size(elek,2)
       % extraiga los nodos del elemento z
       nodz = struct.elements(elek(z),:);

       % busque la posicion del nodo k en el elemento
       posz = find(nodz == k);

       if posz == 1
           % el nodo esta en el nodo 1 local

           % vector tangencial b1
           l1 = n7_n4n(elek(z));
           % vector normal al segmento b1
           b1 = n74cnormal(elek(z),:);

           c17 = sum(u_s7(elek(z),:).*b1).*l1./6;
           c14 = sum(u_s4(elek(z),:).*b1).*l1./4;

           % vector tangencial b2
           l2 = n7_n6n(elek(z));
           % vector normal al segmento b1
           b2 = -n76cnormal(elek(z),:);

           c27 = sum(u_s7(elek(z),:).*b2).*l2./6;
           c26 = sum(u_s6(elek(z),:).*b2).*l2./4;

           vectele = [c17 + c14 + c27 + c26 , c17 + c14 + c27 , c17 + c27 + c26];

       elseif posz == 2
           % el nodo esta en el nodo 2 local

           % vector tangencial b1
           l1 = n7_n4n(elek(z));
           % vector normal al segmento b1
           b1 = -n74cnormal(elek(z),:);

           c17 = sum(u_s7(elek(z),:).*b1).*l1./6;
           c14 = sum(u_s4(elek(z),:).*b1).*l1./4;

           % vector tangencial b2
           l2 = n7_n5n(elek(z));
           % vector normal al segmento b1
           b2 = n75cnormal(elek(z),:);

           c27 = sum(u_s7(elek(z),:).*b2).*l2./6;
           c25 = sum(u_s5(elek(z),:).*b2).*l2./4;

           vectele = [c17 + c14 + c27 , c17 + c14 + c27 + c25 , c27 + c25 + c17];

       elseif posz == 3
           % el nodo esta en el nodo 3 local
                                
           % vector tangencial b1
           l1 = n7_n6n(elek(z));
           % vector normal al segmento b1
           b1 = n76cnormal(elek(z),:);
           
           c17 = sum(u_s7(elek(z),:).*b1).*l1./6;
           c16 = sum(u_s6(elek(z),:).*b1).*l1./4;
                                                       
           % vector tangencial b2
           l2 = n7_n5n(elek(z));
           % vector normal al segmento b1
           b2 = -n75cnormal(elek(z),:);
           
           c27 = sum(u_s7(elek(z),:).*b2).*l2./6;
           c25 = sum(u_s5(elek(z),:).*b2).*l2./4;
           
           vectele = [c17 + c16 + c27 , c27 + c25 + c17 , c17 + c16 + c27 + c25]; 

      end    
       
       matterm2(k,nodz) = vectele + matterm2(k,nodz);
    end
    
end

matterm2 = (matterm2./repmat(dsi,[1 numnodes]));
