% Calcula el gradiente superficial gradsnodes de un campo escalar scfld en
% cada nodo de una malla triangular struct
% basado en Xu Guoliang, "Convergent Discrete Laplace-Beltrami Operators 
% over Triangular Surfaces" 

function gradsnodes = grad_s(struct,scfld)

% preproceso
% ..........
numelements = size(struct.elements,1);
numnodes = size(struct.nodes,1);

% tensor de proyeccion en los elementos
pjele = projtensor(struct.normalele);

    % Calculo del gradiente respecto de Zita1, Zita2 y Zita3 
elementstemp = reshape(struct.elements',numelements*3,1);
scfldtemp = scfld(elementstemp);
scfldtemp = reshape(scfldtemp,3,numelements);
[df1z1, df1z2, df2z1, df2z2, df3z1, df3z2,df1z3] = dfzitalin(scfldtemp);
gradfzita = [df1z1' df1z2' df1z3'];
temp1 = repmat(permute(gradfzita,[3 2 1]),[3 1 1]);
gradtot = sum(struct.jacmat.*temp1,2);
gradtot = repmat(permute(gradtot,[2 1 3]),[3 1 1]);
    % gradiente superficial grad_s(sfld) eN eLeMeNtOS
grads = sum(pjele.*gradtot,2);

% reorganice el vector de gradiente superficial eN eLeMeNtOS.
grads = permute(grads,[1 3 2])';

gradsnodes = zeros(numnodes,3);
for i=1:numelements
    eleind = struct.elements(i,:);
    for k=1:3
        gradsnodes(eleind(k),:) = grads(i,:).*(struct.ds(i)/3) + gradsnodes(eleind(k),:);
    end
end

% gradiente superficial grad_s(gamma) eN LOS NODOS
gradsnodes = gradsnodes./repmat(struct.dsi,[1 3]);
