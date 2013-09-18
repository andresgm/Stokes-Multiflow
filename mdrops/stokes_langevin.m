% CALCULO DEL FLUJO DE STOKES PARA UNA VESICULA
% IMPLEMENTADO FLUJO INFINITO, SEMIINFINITO
% IMPLEMENTADO SINGLE Y DOUBLE LAYER
function [velnode,geom,parms] = stokes_langevin(geom,parms,parmsfem)

numnodes = size(geom.nodes,1);
numelements = size(geom.elements,1);

rkcurv = parms.rkcurv;
rkbend = parms.rkbend;
c0 = parms.c0;
rkgrav = parms.rkgrav;
rkelestat = parms.rkelestat;
rkmaran = parms.rkmaran;

% constante del single layer
rksl = parms.rksl;
% constante del flujo externo
rkextf = parms.rkextf;

gssk = parms.gssk;
csk = parms.csk;

% funcion de green
greenfunction = parms.greenfunction;

% calcule el vector normal a cada nodo
normalandgeoopt.normal = 0;
normalandgeoopt.areas = 1;
normalandgeoopt.vol = 0;
geomprop = normalandgeo(geom,normalandgeoopt);
geom.dsi = geomprop.dsi;
geom.ds = geomprop.ds;
geom.normal = geomprop.normal;
geom.s = geomprop.s;
geom.g = geomprop.g;

polarparms = parms.polarparms;

paropt.tipo = 'extended';
[geom.curv,geom.normal,geom.Kg] = curvparaboloid(geom,paropt);

if rkbend ~= 0
    % calcule el laplace beltrami de la curvatura
    [lapbelmat,geom.Kg] = discretelaplacebeltrami(geom);
    [geom.curv] = curvlb(geom,lapbelmat);
    geom.lapcurv = lapbelmat*geom.curv;
    % geom.lapcurv = lapbel(geom,geom.curv);
    % calcule el delta de fuerza debido a la resistencia al doblamiento.
    % Included spontaneous curvature.
    rdeltafbend = -(rkbend).*...
               ((2.*geom.curv-c0).*...
               (2.*geom.curv.^2-2.*geom.Kg-c0.*geom.curv)+2.*geom.lapcurv);
    geom.rdeltafbend = rdeltafbend;
else
    rdeltafbend = zeros(numnodes,1);
    geom.rdeltafbend = rdeltafbend;
end
    
% Calculo de la tension isotropica Walter et. al.

Mglobal = zeros(3*numnodes,3*numnodes);
rglobal = zeros(3*numnodes,1);

Metricgctr = parmsfem.Metricctr;
Jacmetctr = parmsfem.Jacmetctr;
xg = parmsfem.xg;
wg = parmsfem.wg;

[abase,~,metricg,metricgctr,jacmet,~] = isobase(geom);

tab = ...
    tabskalak(Metricgctr,Jacmetctr,metricg,metricgctr,jacmet,gssk,csk);
chi = zeros(3,3); % chi11, chi22 y chi12 en tres coordenadas cartesianas.

for i = 1:numelements
    for k = 1:7
        Na = p1basis(xg(k,1),xg(k,2));
        dNa = p1dbasis(xg(k,1),xg(k,2));
        wa = wg(k);
        wajac = wa*sqrt(jacmet(i));
        for p = 1:3
            pglob = 3*(geom.elements(i,p)-1)+1;
            chi(1,:) = dNa(p,1).*abase(:,1,i);
            chi(2,:) = dNa(p,2).*abase(:,2,i);
            chi(3,:) = 0.5*(dNa(p,1).*abase(:,2,i)+dNa(p,2).*abase(:,1,i));
            rglobal(pglob) = rglobal(pglob) ...
                + wajac*(chi(1,1)*tab(i,1)+chi(2,1)*tab(i,2)+2*chi(3,1)*tab(i,3));
            rglobal(pglob+1) = rglobal(pglob+1) ...
                + wajac*(chi(1,2)*tab(i,1)+chi(2,2)*tab(i,2)+2*chi(3,2)*tab(i,3));
            rglobal(pglob+2) = rglobal(pglob+2) ...
                + wajac*(chi(1,3)*tab(i,1)+chi(2,3)*tab(i,2)+2*chi(3,3)*tab(i,3));
            Nap = Na(p);
            for q = 1:3
                qglob = 3*(geom.elements(i,q)-1)+1;
                mele = wajac*Nap*Na(q);
                Mglobal(pglob,qglob) = Mglobal(pglob,qglob) + mele;
                Mglobal(pglob+1,qglob+1) = Mglobal(pglob+1,qglob+1) + mele;
                Mglobal(pglob+2,qglob+2) = Mglobal(pglob+2,qglob+2) + mele;
            end
        end
    end
end

qloadv = Mglobal\rglobal;
qload = reshape(qloadv,[3,numnodes])';

% % Esfuerzos normales
% isonorm = repmat(sum(isotens'.*geom.normal,2),[1 3]).*geom.normal;
% % Esfuerzos de marangoni
% rdeltafmaran = isotens'-isonorm;
% 
% 
% %% Resultados
% figure(3);
% grafscfld(geom,normesp(qload));
% hold on
% quiver3(geom.nodes(:,1),geom.nodes(:,2),geom.nodes(:,3),...
%     qload(:,1),qload(:,2),qload(:,3));
% axis equal; view(90,0); xlabel('x1'); ylabel('x2'); zlabel('x3'); colorbar;
% getframe; title('q - skalak');hold off;

rdeltafcurv = rkcurv.*sum(qload.*geom.normal,2);
geom.rdeltafcurv = rdeltafcurv;

% calcule la fuerza de gravedad
if rkgrav ~= 0
    [rdeltafgrav,fuerzagrav] = deltafgrav(geom,rkgrav);
    geom.deltafgrav = rdeltafgrav;
    geom.fuerzagrav = fuerzagrav;
else
    rdeltafgrav = 0;
end

if rkelestat ~= 0
   [rdeltafelestat,fuerzaelest] = deltafelestatic(geom,parms);
   geom.deltafelestat = rdeltafelestat;
   geom.fuerzaelest = fuerzaelest;
else
   rdeltafelestat =0;
end

rdeltaftot = rdeltafcurv + rdeltafgrav + rdeltafbend + rdeltafelestat;

% delta fuerza normal total
geom.rdeltafnorm = rdeltaftot;

% Esfuerzos de marangoni
rdeltafmaran = qload-repmat(sum(qload.*geom.normal,2),[1 3]).*geom.normal;

geom.rdeltafmaran = rkmaran*rdeltafmaran;

% calcule la integral de single layer caso normal y tangencial
rintsln = zeros(numnodes,3);
rintsltanens = zeros(numnodes,3);
rintsltanes = zeros(numnodes,3);

polos = geom.nodes;

for j = 1:numnodes
    
% calcule la integral desingularizada de INT(gij(x,x0)fn(x))ds
    % calcule la funcion de green
    rgreenfcn = stokeslet(polos(j,:),geom.nodes);
    % limpie singularidades
    rgreenfcn(isnan(rgreenfcn)) = 0;
    % calcule deltaf - deltaf*
    rdeltaffpole = repmat((rdeltaftot - rdeltaftot(j)),[1 3]).*geom.normal; 
    % calcule el producto gij*dfi (opt:2)
    rintegrandsl = matvect(rgreenfcn,rdeltaffpole,2);
    % ejecute integral del trapecio
    rintsln(j,:) = inttrapecioa(geom.dsi,rintegrandsl);
    
    % Parte semiinfinita
    if strcmp(parms.flow,'semiinf') == 1
        greenfunctionwall = @greenwall;
        % calcule la funcion de green 
        rgreenwallfcn = stokesletwall(geom.nodes(j,:),geom.nodes,parms.w);
        % limpie singularidades
        rgreenwallfcn(isnan(rgreenwallfcn)) = 0;
        % sume componente pared a stokeslet
        rgreenfcn = rgreenfcn + rgreenwallfcn;
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
        rdeltaffpole = repmat((rdeltaftot - rdeltafim),[1 3]).*geom.normal;
        % calcule el producto gij*dfi (opt:2)
        rintegrandsl = matvect(rgreenwallfcn,rdeltaffpole,2);
        % ejecute integral del trapecio sobre la gota que contiene el polo xj
        rintsln(j,:) = rintsln(j,:) + inttrapecioa(geom.dsi,rintegrandsl);
    end
       
    if rkmaran ~= 0
    % calcule la parte de elementos no singulares de INT(gij(x,x0)ft(x))ds
        % extraiga los elements vecinos al jnodes
        elevecino = geom.element2node{j};
        % calcule la integral no singular
        % calcule el producto gij*ft (opt:2)
        rintegrandsl = matvect(rgreenfcn,rdeltafmaran,2);

        metrica = repmat(geom.g.*(1/3*0.5),[1 3]);
        % integracion de elementos no singulares
        intn1 = rintegrandsl(geom.elements(:,1),:).*metrica;
        intn2 = rintegrandsl(geom.elements(:,2),:).*metrica;
        intn3 = rintegrandsl(geom.elements(:,3),:).*metrica;
        
        intn1(elevecino,:) = 0;
        intn2(elevecino,:) = 0;
        intn3(elevecino,:) = 0;
        rintsltanens(j,:) = sum(intn1 + intn2 + intn3,1);

    % calcule la parte de elementos singulares de INT(gij(x,x0)ft(x))ds
        intval = zeros(1,3);
        for k = elevecino
            % extraiga los indices de los nodos del kelemento
            nodesk = geom.elements(k,:);
            % reorganice los indices del elemento tal que el jnode (polo) quede
            % en P3
            posj = j == nodesk;
            if posj(1) == 1
                % el jnode esta en la posicion 1 local: reorganice
                elenew = [nodesk(2) nodesk(3) nodesk(1)];
            elseif posj(2) == 1
                % el jnode esta en la posicion 2 local: reorganice
                elenew = [nodesk(3) nodesk(1) nodesk(2)];
            end
            triandata.nodes = geom.nodes(elenew,:);
            vectfld = rdeltafmaran(elenew,:);
            triandata.g = geom.g(k);
            intval = intval + ...
                gausslegmatvect2d(greenfunction,vectfld,triandata,polarparms);
            if strcmp(parms.flow,'semiinf') == 1
                intval = intval + ...
              gausslegmatvect2d(greenfunctionwall,vectfld,triandata,polarparms);
            end
        end
        rintsltanes(j,:) = intval; 
    else
            rintsltanens = 0;
            rintsltanes = 0;
    end
        
end

% parte total de INT(gij(x,x0)ft(x))ds
rintslt = rintsltanens + rintsltanes;

% integral total del single layer incluido constante
rintsl = rksl.*(rintslt + rintsln);

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
        [velnode,geom.W] = doublewielandt(geom,rintsl,rextf,geom.velnodeant,geom.W,parms);
    elseif strcmp(parms.dlmod,'subsust') == 1
        % invoque substituciones sucesivas
        [velnode,geom.velnodeant] = sucesivesubs(geom,geom.velnodeant,rintsl,parms);
    end
    geom.velnodeant = velnode;
else
    % calcule la velocidad total
    velnode = rintsl;

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
    if ErrorMax < 1e-6
%         disp(['It Deflaction: ' num2str(zz)]);
%         disp(['Error Deflaction: ' num2str(ErrorMax)]);
        break
    end
end

return
