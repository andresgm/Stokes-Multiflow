% Calcula explicitamente el termino de transporte de un campo escalar
% calcula grad(scfld).vel
function gradscvel = transscfld(struct,scfld,vel)

numelements = size(struct.elements,1);
numnodes = size(struct.nodes,1);

jac = isfield(struct,'jacmat');
dsf = isfield(struct,'ds');
if jac == 0
    % calcule el jacobiano inverso
    ZitaVect = [1/3;1/3];
    [jacomp] = MetricTrans(struct,ZitaVect);
    struct.jacmat = jacomp.jacinv;
end
if dsf == 0
    % calcule las areas
    normalandgeoopt.normal = 1;
    normalandgeoopt.areas = 1;
    normalandgeoopt.vol = 1;
    normalandgeo(struct,normalandgeoopt);
end

% tensor de proyeccion en los elementos
pjele = ProjTensor(struct.normalele);

% Calculo del gradiente respecto de Zita1, Zita2 y Zita3 
ElementsTemp = reshape(struct.elements',numelements*3,1);
scfldtemp = scfld(ElementsTemp);
scfldtemp = reshape(scfldtemp,3,numelements);
[dF1Z1, dF1Z2, dF2Z1, dF2Z2, dF3Z1, dF3Z2,dF1Z3] = dFZitaLin(scfldtemp);

gradFzita = [dF1Z1' dF1Z2' dF1Z3'];

temp1 = repmat(permute(gradFzita,[3 2 1]),[3 1 1]);
gradtot = sum(struct.jacmat.*temp1,2);

gradtot = repmat(permute(gradtot,[2 1 3]),[3 1 1]);
grads = sum(pjele.*gradtot,2);

% reorganice el vector de gradiente superficial.
grads = permute(grads,[1 3 2])';

gradnodes = zeros(numnodes,3);
for i=1:numelements
    eleind = struct.elements(i,:);
    for k=1:3
        gradnodes(eleind(k),:) = grads(i,:).*struct.ds(i) + gradnodes(eleind(k),:);
    end
end

% complete los integrandos del primer termino
gradscvel = (sum(gradnodes.*vel,2)./3)./struct.dsi;

% Calcule el gradiente superficial en los nodos
gs = grad_s(struct,scfld);
% tensor de proyeccion en los nodos
pjnod = ProjTensor(struct.normal);
% refine el gradiente superficial
gs = sum(pjnod.*repmat(permute(gs,[3 2 1]),[3 1 1]),2);
gs = permute(gs,[3 1 2]);
% complete los integrandos del primer termino
gradscvel = sum(gs.*vel,2);
