function StrDotNormal = StrNormDelU(StrFcn,Normal)

NumXpoles = size(StrFcn,1);
StrDotNormal = cell(NumXpoles,1);

[NumElements,NU,NumQuadPoints] = size(Normal);

Temp1 = permute(Normal,[1 3 2]);

    % Cada (NumElements,:) es un la normal en un Punto de Integracion
NormalVect=reshape(Temp1,NumQuadPoints*NumElements,3);

    % Arreglo de valores de la normal para ejecutar producto en paralelo
Temp1 = permute(reshape(NormalVect',1,3*NumQuadPoints*NumElements),[1 3 2]);
NormArray = reshape(repmat(Temp1,[3 3 1]),[3 3 3 NumQuadPoints*NumElements]);
    
for i = 1:NumXpoles
    % Extraiga los tensores de esfuerzo para el iPole
    StrFcnM = StrFcn{i};
    
    % Ejecute Tijk Nk para cada punto de integracion
    StrNormal = sum(StrFcnM.*NormArray,3);
    
   % Permute para que quede array de 3 x 3 x NumQuadPoints*NumElements
    StrNormal = permute(StrNormal,[1 2 4 3]);
    
    StrDotNormal{i} = StrNormal;
end
    
