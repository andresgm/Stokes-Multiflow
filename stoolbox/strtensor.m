function Tst = strtensor(X0,X)

% Numero de Puntos
    NumNodes = size(X,1);

% Asignacion de espacio hipermatriz del tensor de esfuerzos
    Tst = zeros(3,3,3,NumNodes);

% X-X0 para cada componente.
    Xgorro1 = X(:,1) - X0(1);
    Xgorro2 = X(:,2) - X0(2);
    Xgorro3 = X(:,3) - X0(3);

% Vector X-X0
    Xgorro = [Xgorro1 Xgorro2 Xgorro3];
    
% Norma de X-X0 para cada componente y Norma al cubo
    r = NormEsp(Xgorro);
    rquint = r.^5;
    rquint_inv = 1./rquint;

% Quitar Warning de division por cero
    warning off MATLAB:divideByZero

% Calculo de los Xgorro_i*Xgorro_j
    dxx = Xgorro(:,1).^2;
    dxy = Xgorro(:,1).*Xgorro(:,2);
    dxz = Xgorro(:,1).*Xgorro(:,3);
    dyy = Xgorro(:,2).^2;
    dyz = Xgorro(:,2).*Xgorro(:,3);
    dzz = Xgorro(:,3).^2;

% Calculo de las componentes
    cf = (0.75/pi).*rquint_inv;
    
    Tst(1,1,1,:) = dxx.*Xgorro(:,1).*cf;
    Tst(1,1,2,:) = dxy.*Xgorro(:,1).*cf;
    Tst(2,1,1,:) = Tst(1,1,2,:);
    Tst(1,2,1,:) = Tst(1,1,2,:);
  
    Tst(1,1,3,:) = dxz.*Xgorro(:,1).*cf;
    Tst(3,1,1,:) = Tst(1,1,3,:);
    Tst(1,3,1,:) = Tst(1,1,3,:);
    
    Tst(2,1,2,:) = dyy.*Xgorro(:,1).*cf;
    Tst(2,2,1,:) = Tst(2,1,2,:);
    Tst(1,2,2,:) = Tst(2,1,2,:);
    
    Tst(2,1,3,:) = dyz.*Xgorro(:,1).*cf;
    Tst(3,1,2,:) = Tst(2,1,3,:);
    Tst(1,2,3,:) = Tst(2,1,3,:);
    Tst(1,3,2,:) = Tst(2,1,3,:);
    Tst(3,2,1,:) = Tst(2,1,3,:);
    Tst(2,3,1,:) = Tst(2,1,3,:);
   
    Tst(3,1,3,:) = dzz.*Xgorro(:,1).*cf;
    Tst(3,3,1,:) = Tst(3,1,3,:);
    Tst(1,3,3,:) = Tst(3,1,3,:);
    
    Tst(2,2,2,:) = dyy.*Xgorro(:,2).*cf;
    Tst(2,2,3,:) = dyz.*Xgorro(:,2).*cf;
    Tst(3,2,2,:) = Tst(2,2,3,:);
    
    Tst(3,2,3,:) = dzz.*Xgorro(:,2).*cf;
    Tst(3,3,2,:) = Tst(3,2,3,:);
    
    Tst(2,3,2,:) = dyy.*Xgorro(:,3).*cf; 
    Tst(2,3,3,:) = dyz.*Xgorro(:,3).*cf; 
    Tst(3,3,3,:) = dzz.*Xgorro(:,3).*cf;
    
    
