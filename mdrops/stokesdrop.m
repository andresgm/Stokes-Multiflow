% CALCULO DEL FLUJO DE STOKES PARA UNA GOTA
% IMPLEMENTADO FLUJO INFINITO, SEMIINFINITO
% IMPLEMENTADO SINGLE Y DOUBLE LAYER
function [velnode,geom,parms] = stokesdrop(geom,parms)

numnodes = size(geom.nodes,1);
rkgrav = parms.rkgrav;
rkcurv = parms.rkcurv;

rkelect = parms.rkelect;

% constante del single layer
rksl = parms.rksl;
% constante del double layer
rkdl = parms.rkdl;
% constante del flujo externo
rkextf = parms.rkextf;

% calcule el vector normal a cada nodo
normalandgeoopt.normal = 1;
normalandgeoopt.areas = 1;
normalandgeoopt.vol = 1;
geomprop = normalandgeo(geom,normalandgeoopt);

geom.normalele = geomprop.normalele;
geom.normal = geomprop.normal;
geom.dsi = geomprop.dsi;
geom.ds = geomprop.ds;
geom.s = geomprop.s;
geom.vol = geomprop.vol;
geom.jacmat = geomprop.jacmat;

% calculo de la curvatura media
if parms.curvopt == 1
    % calculo mediante paraboloid fitting
    geom.curv = curvparaboloid(geom);
elseif parms.curvopt == 2
    % calculo mediante extended paraboloid fitting/bestparaboloid fitting
    paropt.tipo = 'extended';
    [geom.curv,geom.normal,geom.Kg] = curvparaboloid(geom,paropt);
elseif parms.curvopt == 3
    % calculo mediante laplace - beltrami
    lapbelmat = laplacebeltramimat(geom);
    geom.curv = curvlb(geom,lapbelmat);
end


% calcule la fuerza por tension superficial
if rkcurv ~= 0
    rdeltafcurv = deltafcurv(geom.curv,rkcurv);
else
    rdeltafcurv = 0;
end

% calcule la fuerza de gravedad
if rkgrav ~= 0
    rdeltafgrav = deltafgrav(geom,rkgrav);
    geom.deltafgrav = rdeltafgrav;
else
    rdeltafgrav = 0;
end

% calcule por campo electrico
if rkelect ~= 0
    rdeltafelec = deltafelectric(geom,parms);
    geom.deltafelec = rdeltafelec;
else
    rdeltafelec = 0; 
end

% calcule el delta de fuerza total
rdeltaftot = rdeltafcurv + rdeltafgrav + rdeltafelec;
geom.rdeltafcurv = rdeltafcurv;
geom.rdeltafgrav = rdeltafgrav;
geom.rdeltafelec = rdeltafelec;

% calcule la integral de single layer caso normal y tangencial
rintsln = zeros(numnodes,3);

rdeltafperdrop = cell(geom.numdrops,1);
normalesperdrop = cell(geom.numdrops,1);
dsiperdrop = cell(geom.numdrops,1);
nodesperdrop = cell(geom.numdrops,1);
for k=1:geom.numdrops
    rdeltafperdrop{k} = rdeltaftot(geom.nnodesdrop(k,1):geom.nnodesdrop(k,2),:);  
    normalesperdrop{k} = geom.normal(geom.nnodesdrop(k,1):geom.nnodesdrop(k,2),:);
    dsiperdrop{k} = geom.dsi(geom.nnodesdrop(k,1):geom.nnodesdrop(k,2),:);
    nodesperdrop{k} = geom.nodes(geom.nnodesdrop(k,1):geom.nnodesdrop(k,2),:);
end

for j=1:numnodes
    % calcule la integral de la gota a la que pertenece el polo xj
    % pregunte en que gota esta el polo xj
    for k=1:geom.numdrops
       if j >= geom.nnodesdrop(k,1) && j <= geom.nnodesdrop(k,2)
           gotanum = k;
       else
           gotaext = k; 
       end
    end    
             
    % Calcule la funcion de green
    rgreenfcn = stokeslet(geom.nodes(j,:),nodesperdrop{gotanum});
    % limpie singularidades
    rgreenfcn(isnan(rgreenfcn)) = 0;
    % calcule deltaf - deltaf*
    rdeltaffpole = repmat((rdeltafperdrop{gotanum} - rdeltaftot(j)),[1 3]).*normalesperdrop{gotanum}; 
    % calcule el producto gij*dfi (opt:2)
    rintegrandsl = matvect(rgreenfcn,rdeltaffpole,2);
    % ejecute integral del trapecio sobre la gota que contiene el polo xj
    rintsln(j,:) = inttrapecioa(dsiperdrop{gotanum},rintegrandsl);
  
    % Parte semiinfinita
    if strcmp(parms.flow,'semiinf') == 1
        % calcule la funcion de green 
        rgreenwallfcn = stokesletwall(geom.nodes(j,:),geom.nodes,parms.w);
        % limpie singularidades
        rgreenwallfcn(isnan(rgreenwallfcn)) = 0;
        % determine el deltaf*
        nodex0_im = [geom.nodes(j,1), geom.nodes(j,2), 2*parms.w - geom.nodes(j,3)];
        r = geom.nodes - repmat(nodex0_im,[numnodes 1]);
        rn = normesp(r);
        % Determine el nodo mas cercano a NodeX0_IM
        rdeltafim = rdeltaftot(rn == min(rn));
        if size(rdeltafim,1) ~= 1
            rdeltafim = rdeltafim(1);
        end
        % calcule deltaf - deltaf*
        rdeltaffpole = repmat((rdeltafperdrop{gotanum} - rdeltafim),[1 3]).*normalesperdrop{gotanum};
        % calcule el producto gij*dfi (opt:2)
        rintegrandsl = matvect(rgreenwallfcn,rdeltaffpole,2);
        % ejecute integral del trapecio sobre la gota que contiene el polo xj
        rintsln(j,:) = rintsln(j,:) + inttrapecioa(dsiperdrop{gotanum},rintegrandsl);
    end
    
    % calcule la integral de las demas gotas TODO: ESTA PARA DOS GOTAS - GENERALIZAR PARA J GOTAS
    if geom.numdrops ~= 1
        % Calcule la funcion de green
        [rgreenfcn,closenode] = stokeslet(geom.nodes(j,:),geom.nodes(geom.nnodesdrop(gotaext,1):geom.nnodesdrop(gotaext,2),:));
        % limpie singularidades
        rgreenfcn(isnan(rgreenfcn)) = 0;
        % calcule deltaf - deltaf*
        rdeltafgotaext = rdeltaftot(geom.nnodesdrop(gotaext,1):geom.nnodesdrop(gotaext,2));

        rdeltaffpole = repmat((rdeltafgotaext - rdeltafgotaext(closenode)),[1 3]).*geom.normal(geom.nnodesdrop(gotaext,1):geom.nnodesdrop(gotaext,2),:); 
        % calcule el producto gij*dfi (opt:2)
        rintegrandsl = matvect(rgreenfcn,rdeltaffpole,2);
        % ejecute integral del trapecio de gotaext al polo xj
        rintsln(j,:) = rintsln(j,:) + inttrapecioa(geom.dsi(geom.nnodesdrop(gotaext,1):geom.nnodesdrop(gotaext,2),:),rintegrandsl);
    end
end


rintsln = rksl.*rintsln;

% calcule el flujo externo
if rkextf ~= 0
    rextf = geom.nodes(:,3).*rkextf;
else
    rextf = 0;
end

% calcule wielandt deflaction si lamda .ne. 1
if parms.lamda ~= 1
    if strcmp(parms.dlmod,'deflaction') == 1
        % invoque wielandt deflaction para el double layer
        [velnode,geom.W] = doublewielandt(geom,rintsln,rextf,geom.velnodeant,geom.W,parms);
    elseif strcmp(parms.dlmod,'subsust') == 1
        % invoque substituciones sucesivas
        [velnode,geom.velnodeant] = sucesivesubs(geom,geom.velnodeant,rintsln,parms);
    end
    geom.velnodeant = velnode;
else
    % calcule la velocidad total
    velnode = rintsln;

    % calcule el flujo externo y sumelo
    if rkextf ~= 0
        velnode(:,2) = rextf  + velnode(:,2);
    end
end



return

function [VelAtNode,W] = doublewielandt(geom,IntSingleLayer,UinfAtNode,VelAtNode,W,parms)
% Constante Kappa
Kappa = parms.rkdl*0.5;
    
% se asume siempre que la pared esta en x3 = w = 0
w = 0;
if nargin < 3
% Opciones por defecto
NumSteps = 100;
TolError = 1e-4;
else
    
end
    
Stot = geom.s;
Xpoles = geom.nodes;
dSiV = repmat(geom.dsi,[1 3]);
NumXpoles = size(geom.nodes,1);

% Base Canonica
e = [1 0 0;0 1 0; 0 0 1];
% Surface Centroid Xc
Xc = (1/Stot).*sum(Xpoles.*dSiV,1);
% Xgorro = Xpole - Xc
Xgorro = Xpoles - repmat(Xc,[NumXpoles 1]); 
% Calcule Dij
Xg11 = Xgorro(:,1).^2;
Xg12 = Xgorro(:,1).*Xgorro(:,2);
Xg13 = Xgorro(:,1).*Xgorro(:,3);
Xg22 = Xgorro(:,2).^2;
Xg23 = Xgorro(:,2).*Xgorro(:,3);
Xg33 = Xgorro(:,3).^2;
XgorroQ = normesp(Xgorro).^2;

Dij = zeros(3);
Dij(1,1) = sum((XgorroQ - Xg11).*geom.dsi,1);
Dij(1,2) = sum(-Xg12.*geom.dsi,1);  
Dij(1,3) = sum(-Xg13.*geom.dsi,1);
Dij(2,1) = Dij(1,2);
Dij(2,2) = sum((XgorroQ - Xg22).*geom.dsi,1);
Dij(2,3) = sum(-Xg23.*geom.dsi,1);
Dij(3,1) = Dij(1,3);
Dij(3,2) = Dij(2,3);
Dij(3,3) = sum((XgorroQ - Xg33).*geom.dsi,1);


for zz=1:100
    VelAtNodeAnt = VelAtNode;
    % Calcule IntW = Int(<W(Xpoles) , [ei x Xgorro]>,dS)
    IntW = zeros(3,1);
    for i=1:3
        IntW(i) = sum(sum(W.*cross(repmat(e(i,:),[NumXpoles 1]),Xgorro),2).*geom.dsi,1);
    end
    % Calcule los coeficients Bj = inv(Dij)*IntW
    Bj = Dij\IntW;
    % Calcule Ai = (1/S)*Int(Wi,dS)
    Ai = (1/Stot) .* sum(W.*dSiV,1);
    % W'(Xpole) = Aiei + Biei x Xgorro
    Wprim = repmat(Ai,[NumXpoles,1]) + cross(repmat(Bj',[NumXpoles,1]),Xgorro);

    % Calcule las integrales de double layer.

    % Parte infinita y Total
    IntTstInf = zeros(NumXpoles,3);

    if strcmp(parms.flow,'inf') == 1
        for i=1:NumXpoles
            % Calculo de la integral tensor de esfuerzos parte infinita
            r = geom.nodes - repmat(geom.nodes(i,:),[NumXpoles 1]);
            rn = normesp(r);
            rquint = rn.^5;
            Temp1 = sum(r.*(geom.normal),2);
            Temp2 = W - repmat(W(i,:),[NumXpoles 1]);
            Temp3 = sum(r.*Temp2,2);
            Cf1 = Temp1.*Temp3./rquint;
            Cf1(isnan(Cf1) == 1) = 0;
            % 
            IntTstInf(i,:) = (0.5*W(i,:) + (3/(4*pi)).*sum(dSiV.*repmat(Cf1,[1 3]).*r,1));
        end
            IntTstWall = 0;
    elseif strcmp(parms.flow,'semiinf') == 1
        IntTstWall = zeros(NumXpoles,3);
        for i=1:NumXpoles
             % Extraiga las coordendas del p??lo
            NodeX0 = geom.nodes(i,:);
            NodeX0_IM = [NodeX0(1), NodeX0(2), 2*w - NodeX0(3)];
            R = geom.nodes - repmat(NodeX0_IM,[NumXpoles 1]);
            Rn = normesp(R);
            Rq = Rn.^2;
            Rquint = Rn.^5;
            % Determine el nodo mas cercano a NodeX0_IM
            CloseNode = find(Rn == min(Rn));
            % determine DeltaW = DeltaW - DeltaW0
            if size(CloseNode,1) > 1
                CloseNode = CloseNode(1,1);
            end
            DeltaW = W - repmat(W(CloseNode,:),[NumXpoles 1]);
            RdotN = sum(R.*geom.normal,2);
            RdotW = sum(R.*DeltaW,2);
            WdotN = sum(DeltaW.*geom.normal,2);
            Cf1 = repmat(RdotN.*RdotW./Rquint,[1 3]).*R;
            Part2 = (repmat(5.*(repmat(geom.nodes(i,3),[NumXpoles 1]) - R(:,3)).*(RdotN.*RdotW)./Rq - geom.nodes(i,3).*WdotN,[1 3]).*R...
                + repmat(R(:,3) - repmat(geom.nodes(i,3),[NumXpoles 1]),[1 3]).*(repmat(RdotN, [1 3]).*DeltaW + repmat(RdotW, [1 3]).*geom.normal)...
                + [zeros(NumXpoles,1) zeros(NumXpoles,1) RdotN.*RdotW]);
            Cf2 = -(3*geom.nodes(i,3)/(2*pi))./Rquint;
            Temp1 = Part2.*repmat(Cf2,[1 3]);
            % m=-1 para j=3
            Temp1(:,3) = -Temp1(:,3);
            Part2 = Temp1;
            Part1 = -(3/(4*pi)).*Cf1;
            IntTstWall(i,:) = sum((Part1 + Part2).*dSiV,1);    
        end
    end
    
    IntTstTotal = IntTstInf + IntTstWall;

    % Ensamble de acuerdo a deflaction
    % Calcule la integral Int(<w(Xpole),NormalAtXpole>,dS)
    IntWN = sum(sum(W.*geom.normal,2).*geom.dsi,1);
    % W(Xpole)
    W = parms.rkdl.*(IntTstTotal - Wprim./2 + geom.normal.*(1/(2*Stot)).*IntWN) + IntSingleLayer;
    W(:,2) = W(:,2) + UinfAtNode;
    % Velocidad final de la interfase
    VelAtNode = W + Kappa/(1 - Kappa).*Wprim;

    % Calculo del Error
    ErrorInVel = normesp(VelAtNode - VelAtNodeAnt);
    ErrorMax = max(abs(ErrorInVel));
    if ErrorMax < 1e-3
        disp(['It Deflaction: ' num2str(zz)]);
        disp(['Error Deflaction: ' num2str(ErrorMax)]);
        break
    end
end

return


function [VelAtNode,W] = sucesivesubs(geom,W,IntSingleLayer,parms)

IntSingleLayer = IntSingleLayer./parms.rksl;
KDL = parms.rkdl;
Kgreen = parms.rksl;

 % Calculo de la Integral Cj(X0)
    NumSteps = 100;
    TolError = 1e-3;
    NumXpoles = size(geom.nodes,1);
    dSiV = repmat(geom.dsi,[1 3]);
    VelAtNode = geom.velnodeant;
    for zz=1:NumSteps
        % Calcule las integrales de double layer.

        % Parte infinita y Total
            IntTstInf = zeros(NumXpoles,3);
            IntTstTotal = zeros(NumXpoles,3);
            % W es VelAtNode;
              
        for i=1:NumXpoles
            % Calculo de la integral tensor de esfuerzos parte infinita
            r = geom.nodes - repmat(geom.nodes(i,:),[NumXpoles 1]);
            rn = normesp(r);
            rquint = rn.^5;
            Temp1 = sum(r.*(geom.normal),2);
            Temp2 = W - repmat(W(i,:),[NumXpoles 1]);
            Temp3 = sum(r.*Temp2,2);
            Cf1 = Temp1.*Temp3./rquint;
            Cf1(isnan(Cf1) == 1) = 0;
            % 
            IntTstInf(i,:) = (0.5*W(i,:) + (3/(4*pi)).*sum(dSiV.*repmat(Cf1,[1 3]).*r,1));
            IntTstTotal(i,:) = IntTstInf(i,:);  
        end
        
        % ENSAMBLE LAS INTEGRALES 
        VelAtNodeAnt = VelAtNode;

        VelAtNode = Kgreen.*IntSingleLayer + (KDL).*IntTstTotal;
        
        if parms.rkextf ~= 0
            % Velocidad debido a Uinf
                % Calcule el shear flow en los nodos 
            UinfAtNode = parms.rkextf.*geom.nodes(:,3);
                % Sume el shear flow Uinf pot KGinf
            VelAtNode(:,2) = VelAtNode(:,2) + UinfAtNode;
        end
        W = VelAtNode;
       
        % Calculo del Error
        ErrorInVel = normesp(VelAtNode - VelAtNodeAnt);
        ErrorMax = max(abs(ErrorInVel));
        zz
        if ErrorMax < TolError
            zz
            ErrorMax
            break
        end
            
    end
    
return


% % caso flujo semiinfinito sustituciones sucesivas
% % Extraiga las coordendas del p??lo
% NodeX0 = geom.nodes(i,:);
% NodeX0_IM = [NodeX0(1), NodeX0(2), 2*w - NodeX0(3)];
% R = geom.nodes - repmat(NodeX0_IM,[NumXpoles 1]);
% Rn = NormEsp(R);
% Rq = Rn.^2;
% Rquint = Rn.^5;
% % Determine el nodo mas cercano a NodeX0_IM
% CloseNode = find(Rn == min(Rn));
% % determine DeltaW = DeltaW - DeltaW0
% DeltaW = W - repmat(W(CloseNode,:),[NumXpoles 1]);
% RdotN = sum(R.*geom.normal,2);
% RdotW = sum(R.*DeltaW,2);
% WdotN = sum(DeltaW.*geom.normal,2);
% Cf1 = repmat(RdotN.*RdotW./Rquint,[1 3]).*R;
% Part2 = (repmat(5.*(repmat(geom.nodes(i,3),[NumXpoles 1]) - R(:,3)).*(RdotN.*RdotW)./Rq - geom.nodes(i,3).*WdotN,[1 3]).*R...
%     + repmat(R(:,3) - repmat(geom.nodes(i,3),[NumXpoles 1]),[1 3]).*(repmat(RdotN, [1 3]).*DeltaW + repmat(RdotW, [1 3]).*geom.normal)...
%     + [zeros(NumXpoles,1) zeros(NumXpoles,1) RdotN.*RdotW]);
% Cf2 = (3*geom.nodes(i,3)/(2*pi))./Rquint;
% Temp1 = Part2.*repmat(Cf2,[1 3]);
% Temp1(:,3) = -Temp1(:,3);
% Part2 = Temp1;
% Part1 = -(3/(4*pi)).*(repmat(RdotN.*RdotW./Rquint,[1 3]).*R);
% IntTstTotal(i,:) = IntTstInf(i,:) + sum((Part1 + Part2).*dSiV,1);