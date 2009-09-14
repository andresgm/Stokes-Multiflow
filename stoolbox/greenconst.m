%%Rutina para calcular terminos de la funcion de Green Tij
%%Calcula tambien el nodo mas cercano al polo utilizado
%%x0 = Polo que se esta usando


function [Cf1,r] = greenconst(x0,geom,gotanum,gotaactual,W,pole)

r = geom.nodes(geom.nnodesdrop(gotaactual,1):geom.nnodesdrop(gotaactual,2),:) -...
        repmat(x0,[(geom.nnodesdrop(gotaactual,2)- geom.nnodesdrop(gotaactual,1) + 1) 1]);
rn = normesp(r);
rquint = rn.^5;
% Determina el indice del punto mas cercano a x0
closenode = find(rn == min(rn) & rn > eps);
if size(closenode,1) ~= 1 && size(closenode,1) ~= 0
    closenode = closenode(1);
elseif size(closenode,1) == 0
    % el nodo mas cercano esta en el mismo lugar que x0
    closenode = find(rn == min(rn));
end
Temp1 = sum(r.*(geom.normal(geom.nnodesdrop(gotaactual,1):geom.nnodesdrop(gotaactual,2),:)),2);
if gotanum == gotaactual
    Temp2 = W(geom.nnodesdrop(gotaactual,1):geom.nnodesdrop(gotaactual,2),:) -...
        repmat(W(pole,:),[(geom.nnodesdrop(gotaactual,2)- geom.nnodesdrop(gotaactual,1) + 1) 1]);
else
    Temp2 = W(geom.nnodesdrop(gotaactual,1):geom.nnodesdrop(gotaactual,2),:) -...
        repmat(W(closenode,:),[(geom.nnodesdrop(gotaactual,2)- geom.nnodesdrop(gotaactual,1) + 1) 1]);
end
Temp3 = sum(r.*Temp2,2);
Cf1 = Temp1.*Temp3./rquint;
Cf1(isnan(Cf1) == 1) = 0;