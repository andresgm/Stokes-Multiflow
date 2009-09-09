function IntDoubleLayer = StrNormDelU2(StrFcn,Normal,VelAtIntPtV,VelAtNode,IndexTable,NumElements,RedInv)

NumXpoles = size(StrFcn,1);
IntDoubleLayer = cell(NumXpoles,1);
CellTemp = cell(1,1);

[m,NU,o] = size(Normal);
NumQuadPoints = size(VelAtIntPtV,1)/NumElements;

if o > 1  
    % Entrada tipo array Dim(m,3,o)
    Temp1 = permute(Normal,[1 3 2]);
        % Cada (NumElements,:) es un la normal en un Punto de Integracion
    NormalVect=reshape(Temp1,NumQuadPoints*NumElements,3);
    % Arreglo de valores de la normal para ejecutar producto en paralelo
    Temp1 = permute(reshape(NormalVect',1,3*NumQuadPoints*NumElements),[1 3 2]);
    NormArray = reshape(repmat(Temp1,[3 3 1]),[3 3 3 NumQuadPoints*NumElements]);
else
    % Entrada tipo lista
    NormalVect = Normal;
    % Arreglo de valores de la normal para ejecutar producto en paralelo
    Temp1 = permute(reshape(NormalVect',1,3*m),[1 3 2]);
    NormArray = reshape(repmat(Temp1,[3 3 1]),[3 3 3 m]);
end

    
    
for i = 1:NumXpoles
    % Extraiga los tensores de esfuerzo para el iPole
    StrFcnM = StrFcn{i};
    % Extraiga la velocidad del iPole
    VelPole = VelAtNode(IndexTable(i),:);
    % Realice U(X) - U(X0)
    DeltaVel = VelAtIntPtV - repmat(VelPole,[NumElements*NumQuadPoints,1]);
    % Reconstruya DeltaVel a dim(NumElements,3,NumQuadPoints)
%     Temp1=reshape(DeltaVel,NumElements,NumQuadPoints,3);
%     DeltaVel = permute(Temp1,[1 3 2]);
    if size(RedInv,1) > 1
        % Contraiga el array DeltaVel
        DeltaVel = DeltaVel(RedInv,:);
    end
    % Ejecute Tijk Nk para cada punto de integracion
    StrNormal = sum(StrFcnM.*NormArray,3);
    
   % Permute para que quede array de 3 x 3 x NumQuadPoints*NumElements
    StrNormal = permute(StrNormal,[1 2 4 3]);
    CellTemp{1} = StrNormal;
    
   % Invoque GreenDotDeltaF para ejecutar Dij*Vi
%     IntDoubleLayerT = GreenDotDeltaF(CellTemp,DeltaVel);
    IntDoubleLayerT = matvect(StrNormal,DeltaVel,2);
    
    IntDoubleLayer{i} = IntDoubleLayerT;
end
    
