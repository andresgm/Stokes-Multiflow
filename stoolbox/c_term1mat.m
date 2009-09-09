% calcula la matriz del termino 1 de sulfactantes

function matterm1 = c_term1mat(struct,velt)

numnodes = size(struct.nodes,1);
numelements = size(struct.elements,1);

% OJO CON ESTO BORRAR... SOLO PARA PRUEBA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% w = zeros(numnodes,3);
% vel = zeros(numnodes,3);
% velt = zeros(numnodes,3);;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(struct,'normalele') ~= 1 || isfield(struct,'normal') ~= 1
   [normalnode,normalele] = normal(struct); 
else
    normalele = struct.normalele;
    normalnode = struct.normal;
end

if isfield(struct,'jacmat') ~= 1
    % no tiene el inverso del jacobiano calculelo
    jacomp = metrictrans(struct,[1/3;1/3]);
    invjac = jacomp.jacinv;
else
    invjac = struct.jacmat;
end

if isfield(struct,'ds') ~= 1
    % no hay areas de cada elemento calculadas
    [s,ds,dsi] = areas(struct);
else
    dsi = struct.dsi;
    ds = struct.ds;
end

% pjnod = ProjTensor(struct.normal);
pjele = projtensor(normalele);

% Velocidad tangencial
velt = velt - repmat(sum(velt.*normalnode,2),[1 3]).*normalnode;

rk = zeros(3,3,numelements);
parfor k=1:numelements
    rk(:,:,k) = pjele(:,:,k)*invjac(:,:,k).*(ds(k)/3);
end

% matriz gradiente superficial en elementos
rnk = rk;
rnk(:,3,:) = -rk(:,1,:)-rk(:,2,:)+rk(:,3,:);

% armar la matriz global en nodos
% tres matrices para cada componente

mat1 = zeros(numnodes,numnodes);
mat2 = zeros(numnodes,numnodes);
mat3 = zeros(numnodes,numnodes);
for k=1:numelements
    % extraiga los elementos vecinos al jnode
    nodoele = struct.elements(k,:);
    for h=1:3
        mat1(nodoele(h),nodoele(1)) = rnk(1,1,k) + mat1(nodoele(h),nodoele(1));
        mat1(nodoele(h),nodoele(2)) = rnk(1,2,k) + mat1(nodoele(h),nodoele(2));
        mat1(nodoele(h),nodoele(3)) = rnk(1,3,k) + mat1(nodoele(h),nodoele(3));
        
        mat2(nodoele(h),nodoele(1)) = rnk(2,1,k) + mat2(nodoele(h),nodoele(1));
        mat2(nodoele(h),nodoele(2)) = rnk(2,2,k) + mat2(nodoele(h),nodoele(2));
        mat2(nodoele(h),nodoele(3)) = rnk(2,3,k) + mat2(nodoele(h),nodoele(3));
        
        mat3(nodoele(h),nodoele(1)) = rnk(3,1,k) + mat3(nodoele(h),nodoele(1));
        mat3(nodoele(h),nodoele(2)) = rnk(3,2,k) + mat3(nodoele(h),nodoele(2));
        mat3(nodoele(h),nodoele(3)) = rnk(3,3,k) + mat3(nodoele(h),nodoele(3));
    end
end

t1 = mat1.*repmat(velt(:,1),[1 numnodes]);
t2 = mat2.*repmat(velt(:,2),[1 numnodes]);
t3 = mat3.*repmat(velt(:,3),[1 numnodes]);
matterm1 = sparse((t1 + t2 + t3)./repmat(dsi,[1 numnodes]));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gradsmat = zeros(numelements,3);
% for k=1:numelements
%     % extraiga los nodos de cada elemento
%     nodoele = geom.elements(k,:);
%     % arme la derivada df/x1;
%     duz = [gamma(nodoele(1)) - gamma(nodoele(3)); ...
%         gamma(nodoele(2)) - gamma(nodoele(3)); ...
%         gamma(nodoele(3))];
%     gradsmat(k,:) = (rk(:,:,k)*duz)';
% end
% gradsnodes2 = zeros(numnodes,3);
% for i=1:numelements
%     eleind = geom.elements(i,:);
%     for k=1:3
%         gradsnodes2(eleind(k),:) = gradsmat(i,:) + gradsnodes2(eleind(k),:);
%     end
% end
% 
% % gradiente superficial grad_s(gamma) En LOS nODOS
% gradsnodes2 = gradsnodes2./repmat(geom.dsi,[1 3]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gradsmat3 = zeros(numelements,3);
% for k=1:numelements
%     % extraiga los nodos de cada elemento
%     nodoele = geom.elements(k,:);
%     gradsmat3(k,:) = (rnk(:,:,k)*gamma(nodoele))';
% end
% 
% gradsnodes3 = zeros(numnodes,3);
% for i=1:numelements
%     eleind = geom.elements(i,:);
%     for k=1:3
%         gradsnodes3(eleind(k),:) = gradsmat3(i,:) + gradsnodes3(eleind(k),:);
%     end
% end
% 
% % gradiente superficial grad_s(gamma) En LOS nODOS
% gradsnodes3 = gradsnodes3./repmat(geom.dsi,[1 3]);
