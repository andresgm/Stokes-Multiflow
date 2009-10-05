% CALCULO DEL FLUJO DE STOKES PARA N GOTAS
% IMPLEMENTADO FLUJO INFINITO
% IMPLEMENTADO SINGLE LAYER
function [velnode,geom] = stokesmdrop(geom,parms)

numnodes = size(geom.nodes,1);
rkgrav = parms.rkgrav;
rkcurv = parms.rkcurv;

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

% calcule la curvatura media
% [geom.curv,geom.normal] = curvbestparaboloid(geom);
geom.curv = curvparaboloid(geom);
% geom.curv = curvlb(geom);

% calcule los esfuerzos por tension superficial
rdeltafcurv = deltafcurv(geom.curv,rkcurv);

% calcule la fuerza de gravedad
if rkgrav ~= 0
    rdeltafgrav = deltafgrav(geom.nodes,rkgrav);
else
    rdeltafgrav = 0;
end

% calcule el delta de fuerza total NORMAL
rdeltaftot = rdeltafcurv + rdeltafgrav;

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

nnodesdrop = geom.nnodesdrop;
numdrops = geom.numdrops;
nodes = geom.nodes;
normales = geom.normal;
dsi = geom.dsi;
for j=1:numnodes
    % extraiga el polo xj
    % calcule la integral de la gota a la que pertenece el polo xj
        % pregunte en que gota esta el polo xj
        for k=1:numdrops
           if j >= nnodesdrop(k,1) && j <= nnodesdrop(k,2)
               gotanum = k;
           end
        end
    
    % Calcule la funcion de green
    rgreenfcn = stokeslet(nodes(j,:),nodesperdrop{gotanum});
    % limpie singularidades
    rgreenfcn(isnan(rgreenfcn)) = 0;
    % calcule deltaf - deltaf*
    rdeltaffpole = repmat((rdeltafperdrop{gotanum} - rdeltaftot(j)),[1 3]).*normalesperdrop{gotanum}; 
    % calcule el producto gij*dfi (opt:2)
    rintegrandsl = matvect(rgreenfcn,rdeltaffpole,2);
    % ejecute integral del trapecio sobre la gota que contiene el polo xj
    rintsln(j,:) = inttrapecioa(dsiperdrop{gotanum},rintegrandsl);
    
    
    % calcule la integral de las demas gotas
    if numdrops ~= 1
        for m=1:geom.numdrops
            if m ~= gotanum
                % Calcule la funcion de green 
                [rgreenfcn,closenode] = stokeslet(nodes(j,:),nodes(nnodesdrop(m,1):nnodesdrop(m,2),:));
                % limpie singularidades
                rgreenfcn(isnan(rgreenfcn)) = 0;
                % calcule deltaf - deltaf*
                rdeltafgotaext = rdeltaftot(nnodesdrop(m,1):nnodesdrop(m,2));
                rdeltaffpole = repmat((rdeltafgotaext - rdeltafgotaext(closenode)),[1 3]).*normales(nnodesdrop(m,1):nnodesdrop(m,2),:); 
                % calcule el producto gij*dfi (opt:2)
                rintegrandsl = matvect(rgreenfcn,rdeltaffpole,2);
                % ejecute integral del trapecio de gotaext al polo xj
                rintsln(j,:) = rintsln(j,:) + inttrapecioa(dsi(nnodesdrop(m,1):nnodesdrop(m,2),:),rintegrandsl);
            else
            end
        end
    end
end



rintsln = rksl.*rintsln;

% % calcule el flujo externo
% rextf = geom.nodes(:,3).*rkextf;
% 
% % calcule la velocidad total
% velnode = rintsln;
% velnode(:,2) = rextf  + velnode(:,2);

% calcule el flujo externo
if rkextf ~= 0
    rextf = geom.nodes(:,3).*rkextf;
else
    rextf = 0;
end

% calcule wielandt deflaction si lamda .ne. 1
if parms.lamda ~= 1
    % invoque wielandt deflaction para el double layer
        [velnode,geom.W] = doublewielandt(geom,rintsln,rextf,geom.velnodeant,geom.W,parms);
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
Xc = zeros(geom.numdrops,3);
Xgorro = cell(geom.numdrops,1);
Dij = cell(geom.numdrops,1);

% Base Canonica
e = [1 0 0;0 1 0; 0 0 1];
% Surface Centroid Xc

for j=1:geom.numdrops
Xc(j,:) = (1/Stot(j)).*sum(Xpoles(geom.nnodesdrop(j,1):geom.nnodesdrop(j,2),:).*dSiV(geom.nnodesdrop(j,1):geom.nnodesdrop(j,2),:),1);
% Xgorro = Xpole - Xc
Xgorro{j} = Xpoles(geom.nnodesdrop(j,1):geom.nnodesdrop(j,2),:) - repmat(Xc(j,:),[(geom.nnodesdrop(j,2)- geom.nnodesdrop(j,1) + 1) 1]);

% Calcule Dij
Xg11 = Xgorro{j}(:,1).^2;
Xg12 = Xgorro{j}(:,1).*Xgorro{j}(:,2);
Xg13 = Xgorro{j}(:,1).*Xgorro{j}(:,3);
Xg22 = Xgorro{j}(:,2).^2;
Xg23 = Xgorro{j}(:,2).*Xgorro{j}(:,3);
Xg33 = Xgorro{j}(:,3).^2;
XgorroQ = normesp(Xgorro{j}).^2;

Dij{j}(1,1) = sum((XgorroQ - Xg11).*geom.dsi(geom.nnodesdrop(j,1):geom.nnodesdrop(j,2),:),1);
Dij{j}(1,2) = sum(-Xg12.*geom.dsi(geom.nnodesdrop(j,1):geom.nnodesdrop(j,2),:),1);  
Dij{j}(1,3) = sum(-Xg13.*geom.dsi(geom.nnodesdrop(j,1):geom.nnodesdrop(j,2),:),1);
Dij{j}(2,1) = Dij{j}(1,2);
Dij{j}(2,2) = sum((XgorroQ - Xg22).*geom.dsi(geom.nnodesdrop(j,1):geom.nnodesdrop(j,2),:),1);
Dij{j}(2,3) = sum(-Xg23.*geom.dsi(geom.nnodesdrop(j,1):geom.nnodesdrop(j,2),:),1);
Dij{j}(3,1) = Dij{j}(1,3);
Dij{j}(3,2) = Dij{j}(2,3);
Dij{j}(3,3) = sum((XgorroQ - Xg33).*geom.dsi(geom.nnodesdrop(j,1):geom.nnodesdrop(j,2),:),1);
end

for zz=1:1000
    VelAtNodeAnt = VelAtNode;
    % Calcule IntW = Int(<W(Xpoles) , [ei x Xgorro]>,dS)
    for j=1:geom.numdrops
        IntW = zeros(3,1);
        for i=1:3
            IntW(i) = sum(sum(W(geom.nnodesdrop(j,1):geom.nnodesdrop(j,2),:).*cross(repmat(e(i,:),[(geom.nnodesdrop(j,2)- ...
                geom.nnodesdrop(j,1) + 1) 1]),Xgorro{j}),2).*geom.dsi(geom.nnodesdrop(j,1):geom.nnodesdrop(j,2),:),1);
        end
        % Calcule los coeficients Bj = inv(Dij)*IntW
        Bj(:,j) = Dij{j}\IntW;
        % Calcule Ai = (1/S)*Int(Wi,dS)
        Ai(j,:) = (1/Stot(j)) .* sum(W(geom.nnodesdrop(j,1):geom.nnodesdrop(j,2),:).*dSiV(geom.nnodesdrop(j,1):geom.nnodesdrop(j,2),:),1);
        % W'(Xpole) = Aiei + Biei x Xgorro
        Wprim(geom.nnodesdrop(j,1):geom.nnodesdrop(j,2),:) = repmat(Ai(j,:),[(geom.nnodesdrop(j,2)- ...
            geom.nnodesdrop(j,1) + 1),1]) + cross(repmat(Bj(:,j)',[(geom.nnodesdrop(j,2) - geom.nnodesdrop(j,1) + 1),1]),Xgorro{j});
    end
    % Calcule las integrales de double layer.

    % Parte infinita y Total
    IntTstInf = zeros(NumXpoles,3);
    for i=1:NumXpoles
        % Pregunte en que gota esta el polo i
        for k=1:geom.numdrops
           if i >= geom.nnodesdrop(k,1) && i <= geom.nnodesdrop(k,2)
               gotanum = k;
           end
        end
        for j=1:geom.numdrops
                [Cf1,r] = greenconst(geom.nodes(i,:),geom,gotanum,j,W,i);
                if j==gotanum
                    IntTstInf(i,:) = IntTstInf(i,:) + (0.5*W(i,:) + (3/(4*pi)).*...
                        sum((dSiV(geom.nnodesdrop(j,1):geom.nnodesdrop(j,2),:).*repmat(Cf1,[1 3])).*r,1));
                else
                    IntTstInf(i,:) = IntTstInf(i,:) + ((3/(4*pi)).*...
                    sum((dSiV(geom.nnodesdrop(j,1):geom.nnodesdrop(j,2),:).*repmat(Cf1,[1 3])).*r,1));
                end
        end
    end


    % Ensamble de acuerdo a deflaction
    % Calcule la integral Int(<w(Xpole),NormalAtXpole>,dS)
    for j=1:geom.numdrops
        IntWN(j) = sum(sum(W(geom.nnodesdrop(j,1):geom.nnodesdrop(j,2),:).*...
        geom.normal(geom.nnodesdrop(j,1):geom.nnodesdrop(j,2),:),2).*geom.dsi(geom.nnodesdrop(j,1):geom.nnodesdrop(j,2),:),1);
        norsubS(geom.nnodesdrop(j,1):geom.nnodesdrop(j,2),:) = geom.normal(geom.nnodesdrop(j,1):geom.nnodesdrop(j,2),:).*(1/(2*Stot(j)));
        % W(Xpole)
        W(geom.nnodesdrop(j,1):geom.nnodesdrop(j,2),:) = parms.rkdl.*(IntTstInf(geom.nnodesdrop(j,1):geom.nnodesdrop(j,2),:)...
            - Wprim(geom.nnodesdrop(j,1):geom.nnodesdrop(j,2),:)./2 + norsubS(geom.nnodesdrop(j,1):geom.nnodesdrop(j,2),:).*IntWN(j))...
            + IntSingleLayer(geom.nnodesdrop(j,1):geom.nnodesdrop(j,2),:);
    end
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
