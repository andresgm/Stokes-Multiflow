function VelAtNode = integrationf(data,coordinates)
      

%% Calculo de la integral de Single Layer

numpoles = size(coordinates,1);

% calcule los esfuerzos por tension superficial
rdeltafcurv = deltafcurv(data.geom.curv,data.parms.rkcurv);

% calcule la fuerza de gravedad
if data.parms.rkgrav ~= 0
    rdeltafgrav = deltafgrav(data.geom.nodes,data.parms.rkgrav);
else
    rdeltafgrav = 0;
end

% calcule el delta de fuerza total NORMAL
rdeltaftot = rdeltafcurv + rdeltafgrav;

rintsln = zeros(numpoles,3);

rdeltafperdrop = cell(data.geom.numdrops,1);
normalesperdrop = cell(data.geom.numdrops,1);
dsiperdrop = cell(data.geom.numdrops,1);
nodesperdrop = cell(data.geom.numdrops,1);
for k=1:data.geom.numdrops
rdeltafperdrop{k} = rdeltaftot(data.geom.nnodesdrop(k,1):data.geom.nnodesdrop(k,2),:);  
normalesperdrop{k} = data.geom.normal(data.geom.nnodesdrop(k,1):data.geom.nnodesdrop(k,2),:);
dsiperdrop{k} = data.geom.dsi(data.geom.nnodesdrop(k,1):data.geom.nnodesdrop(k,2),:);
nodesperdrop{k} = data.geom.nodes(data.geom.nnodesdrop(k,1):data.geom.nnodesdrop(k,2),:);
end

for j=1:numpoles
    % extraiga el polo xj
    % calcule la integral de la gota a la que pertenece el polo xj
        % pregunte en que gota esta el polo xj
        for k=1:data.geom.numdrops
           if j >= data.geom.nnodesdrop(k,1) && j <= data.geom.nnodesdrop(k,2)
               gotanum = k;
           end
        end
    
    % Calcule la funcion de green
    rgreenfcn = stokeslet(coordinates(j,:),nodesperdrop{gotanum});
    % limpie singularidades
    rgreenfcn(isnan(rgreenfcn)) = 0;
    % calcule deltaf - deltaf*
    rdeltaffpole = repmat(rdeltafperdrop{gotanum},[1 3]).*normalesperdrop{gotanum}; 
    % calcule el producto gij*dfi (opt:2)
    rintegrandsl = matvect(rgreenfcn,rdeltaffpole,2);
    % ejecute integral del trapecio sobre la gota que contiene el polo xj
    rintsln(j,:) = inttrapecioa(dsiperdrop{gotanum},rintegrandsl);
    
    
    % calcule la integral de las demas gotas
    if data.geom.numdrops ~= 1
        for m=1:data.geom.numdrops
            if m ~= gotanum
                % Calcule la funcion de green 
                [rgreenfcn] = stokeslet(coordinates(j,:),data.geom.nodes(data.geom.nnodesdrop(m,1):data.geom.nnodesdrop(m,2),:));
                % limpie singularidades
                rgreenfcn(isnan(rgreenfcn)) = 0;
                % calcule deltaf - deltaf*
                rdeltafgotaext = rdeltaftot(data.geom.nnodesdrop(m,1):data.geom.nnodesdrop(m,2));
                rdeltaffpole = repmat(rdeltafgotaext,[1 3]).*data.geom.normal(data.geom.nnodesdrop(m,1):data.geom.nnodesdrop(m,2),:); 
                % calcule el producto gij*dfi (opt:2)
                rintegrandsl = matvect(rgreenfcn,rdeltaffpole,2);
                % ejecute integral del trapecio de gotaext al polo xj
                rintsln(j,:) = rintsln(j,:) + inttrapecioa(data.geom.dsi(data.geom.nnodesdrop(m,1):data.geom.nnodesdrop(m,2),:),rintegrandsl);
            else
            end
        end
    end
end

% calcule el flujo externo
if data.parms.rkextf ~= 0
    rextf = coordinates(:,3).*data.parms.rkextf;
else
    rextf = 0;
end

% calcule wielandt deflaction si lamda .ne. 1
if data.parms.lamda ~= 1
    % invoque wielandt deflaction para el double layer
        [VelAtNode] = doublelayer(data.geom,rintsln,rextf,data.velnode,data.parms,coordinates);
else
    % calcule la velocidad total
    VelAtNode = rintsln;
    % calcule el flujo externo y sumelo
    if data.parms.rkextf ~= 0
        VelAtNode(:,2) = rextf  + VelAtNode(:,2);
    end
end

return

%% Calculo de la integral de Double Layer

function [VelAtNode] = doublelayer(geom,IntSingleLayer,UinfAtNode,velnode,parms,coordinates)
% Constante Kappa
 Kappa = parms.rkdl*0.5*(parms.lamda+1);
dSiV = repmat(geom.dsi,[1 3]);

    IntTstInf = zeros(size(coordinates,1),3);
    for i=1:size(coordinates,1)
        for j=1:geom.numdrops
            r = geom.nodes(geom.nnodesdrop(j,1):geom.nnodesdrop(j,2),:) -...
                repmat(coordinates(i,:),[(geom.nnodesdrop(j,2)- geom.nnodesdrop(j,1) + 1) 1]);
            rn = normesp(r);
            rquint = rn.^5;
            Temp1 = sum(r.*(geom.normal(geom.nnodesdrop(j,1):geom.nnodesdrop(j,2),:)),2);
            Temp2 = velnode(geom.nnodesdrop(j,1):geom.nnodesdrop(j,2),:);
            Temp3 = sum(r.*Temp2,2);
            Cf1 = Temp1.*Temp3./rquint;
            Cf1(isnan(Cf1) == 1) = 0;
            IntTstInf(i,:) = IntTstInf(i,:) + ((3/(4*pi)).*...
                sum((dSiV(geom.nnodesdrop(j,1):geom.nnodesdrop(j,2),:).*repmat(Cf1,[1 3])).*r,1));
        end
    end

    VelAtNode = Kappa.*IntTstInf + IntSingleLayer;
    VelAtNode(:,2) = VelAtNode(:,2) + UinfAtNode;
    
return
