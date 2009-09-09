% CALCULO DEL FLUJO DE STOKES PARA 2 GOTAS
% IMPLEMENTADO FLUJO INFINITO
% IMPLEMENTADO SINGLE LAYER
function [velnode,geom] = stokes(geom,parms)

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
parfor j=1:numnodes
    % extraiga el polo xj
    
    % calcule la integral de la gota a la que pertenece el polo xj
        % pregunte en que gota esta el polo xj
        for k=1:numdrops
           if j >= nnodesdrop(k,1) && j <= nnodesdrop(k,2)
               gotanum = k;
           else
               gotaext = k; 
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
    
    
    % calcule la integral de las demas gotas TODO: ESTA PARA DOS GOTAS - GENERALIZAR PARA J GOTAS
    if numdrops ~= 1
        % Calcule la funcion de green
        [rgreenfcn,closenode] = stokeslet(nodes(j,:),nodes(nnodesdrop(gotaext,1):nnodesdrop(gotaext,2),:));
        % limpie singularidades
        rgreenfcn(isnan(rgreenfcn)) = 0;
        % calcule deltaf - deltaf*
        rdeltafgotaext = rdeltaftot(nnodesdrop(gotaext,1):nnodesdrop(gotaext,2));

        rdeltaffpole = repmat((rdeltafgotaext - rdeltafgotaext(closenode)),[1 3]).*normales(nnodesdrop(gotaext,1):nnodesdrop(gotaext,2),:); 
        % calcule el producto gij*dfi (opt:2)
        rintegrandsl = matvect(rgreenfcn,rdeltaffpole,2);
        % ejecute integral del trapecio de gotaext al polo xj
        rintsln(j,:) = rintsln(j,:) + inttrapecioa(dsi(nnodesdrop(gotaext,1):nnodesdrop(gotaext,2),:),rintegrandsl);
    end
end



rintsln = rksl.*rintsln;

% calcule el flujo externo
rextf = geom.nodes(:,3).*rkextf;

% calcule la velocidad total
velnode = rintsln;
velnode(:,2) = rextf  + velnode(:,2);

