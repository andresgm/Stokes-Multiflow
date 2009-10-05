% RUTINA DE PARABOLOID FITTING DADA POR:
% Curvature Estimation for unstructured triangulations of Surfaces
% R.V Garimella & B.K Swartz 
% technical Report los alamos national laboratory
% equation z = ax^2 + bxY + cY^2
% Calcula la curvatura media y curvatura gaussiana

function  [cmean,normal_f,Kg] = curvparaboloid(geom,paropt)

if nargin < 2
    tipo = 'single';
else
    tipo = paropt.tipo;
    if strcmp(tipo,'extended') == 1 && isfield(paropt,'tol') ~= 1 ...
           && isfield(paropt,'itmax') ~= 1
        % use opciones para extended por defecto

        tol = 1e-3;
        itmax = 10;
    elseif strcmp(tipo,'extended') == 1 && isfield(paropt,'tol') == 1 ...
           && isfield(paropt,'itmax') == 1
        tol = paropt.tol;
        itmax = paropt.itmax;
    end       
end

% recupere contadores.
numnodes = size(geom.nodes,1);

if isfield(geom,'nodecon2node') ~= 1
    % calcule la conectividad de nodos
    nodecon2node = node2node(geom.elements);   
else
    nodecon2node = geom.nodecon2node;    
end

if  isfield(geom,'normal') ~= 1
    % no se ha calculado la normal inicial calcularla
    norm_i = normal(geom);
else
    norm_i = geom.normal;
end       

if strcmp(tipo,'single') == 1
% calcule el tensor de proyeccion tangente a los nodos
pg = projtensor(norm_i);    
end

% defina espacio en memoria
cmean = zeros(numnodes,1);
Kg = zeros(numnodes,1);
normalbp = zeros (numnodes,3);

if strcmp(tipo,'single') == 1
    % r1 es extraer la primera columna del tensor de proyeccion
    r1 = permute(pg(:,1,:),[3 1 2]);
    r1 = univect(r1);
    r2 = crossv(norm_i,r1);
end

for i=1:numnodes
    % extraiga los indices de los nodos vertices conectados al i
    nodesadj = nodecon2node{i};

    % calcule las posiciones locales
    % extraiga las coordenadas de los nodos conectados al is
    xglobnodes = geom.nodes(nodesadj,:)';
    nodecen = geom.nodes(i,:)';

    deltax = zeros(3,size(nodesadj,1));
    for k=1:size(nodesadj,1)
        temp1 = xglobnodes(:,k) - nodecen;
        deltax(:,k) = temp1;
    end

    if strcmp(tipo,'single') == 1
        % paraboloid fitting based in z = ax^2 + bxY + cY^2
        mr = [r1(i,:);r2(i,:);norm_i(i,:)];
        xlocnodes = zeros(3,size(nodesadj,1));
        for k=1:size(nodesadj,1)
            xlocnodes(:,k) = mr*deltax(:,k);
        end
        xlocnodes = xlocnodes';
        [cmean(i),Kg(i)] = simplequadric(xlocnodes);
        
    elseif strcmp(tipo,'extended') == 1
         % paraboloid fitting based in z = ax^2 + bxY + cY^2 + dx + eY
        [cmean(i),normalbp(i,:),Kg(i)] = extendedquatric(deltax,norm_i(i,:),tol,itmax);
    end

end

if strcmp(tipo,'single') == 1
    % Si se usa metodo 'single' la normal no se actualiza
    normal_f = normal_i;
elseif strcmp(tipo,'extended') == 1
    normal_f = normalbp;
end

end
   
function [cmeancen,Kg] = simplequadric(xlocnodes)
       
% ensamble la matriz del sistema lineal
    col1 = xlocnodes(:,1).^2;
    col2 = xlocnodes(:,1).*xlocnodes(:,2);
    col3 = xlocnodes(:,2).^2;
    matriz = [col1 col2 col3];
    zeta = xlocnodes(:,3);
    
    % Solucion en minimos cuadrados
    alfa = matriz\zeta;

    % curvatura media
    cmeancen = -(alfa(1) + alfa(3));
    
%     % curvaturas principales
%     k1 = -(alfa(1) + alfa(3) + sqrt((alfa(1) - alfa(3))^2 + alfa(2)^2));
%     k2 = -(alfa(1) + alfa(3) - sqrt((alfa(1) - alfa(3))^2 + alfa(2)^2));
    
    % curvatura gaussiana
    Kg = 4.*alfa(1).*alfa(3) - alfa(2)^2;
    
end

function [cmeancen,normalcen,Kg] = extendedquatric(deltax,normalnew,tol,itmax)

    normalant = [0,0,0];
    it = 0;
    while normesp(normalnew - normalant) > tol
        it = it + 1;
        % matriz de rotaciones
        r3 = normalnew;
        pga = projtensor(normalnew);
        r1 = pga(:,1)';
        r1 = univect(r1);
        r2 = crossv(r3,r1);
            
        mr = [r1;r2;r3];
        
        xlocnodes = zeros(3,size(deltax,2));
        for k=1:size(deltax,2)
            xlocnodes(:,k) = mr*deltax(:,k);
        end
        xlocnodes = xlocnodes';
        
    % ensamble la matriz del sistema lineal
        col1 = xlocnodes(:,1).^2;
        col2 = xlocnodes(:,1).*xlocnodes(:,2);
        col3 = xlocnodes(:,2).^2;
        col4 = xlocnodes(:,1);
        col5 = xlocnodes(:,2);
        
        matriz = [col1 col2 col3 col4 col5];
        zeta = xlocnodes(:,3);

        % Solucion en minimos cuadrados
        alfa = matriz\zeta;

        % extraccion de propiedades
        % normal en la base local
        normalant = normalnew;
        normalnewl = (1/(alfa(4)^2 + alfa(5)^2 + 1)^0.5).*[-alfa(4),-alfa(5),1]';
        % transforme la normal a la base Global
        normalnew = (mr'*normalnewl)';
        % curvatura media
        cmeancen = -(alfa(1) + alfa(3) + alfa(1)*alfa(5)^2 + alfa(3)*alfa(4)^2 ...
            - alfa(2)*alfa(4)*alfa(5))/(1 + alfa(4)^2 + alfa(5)^2)^1.5;
        
        
        % curvatura gaussiana
        Kg= ((4.*alfa(1).*alfa(3) - alfa(2)^2)./(1 + alfa(4)^2 + alfa(5)^2)^2);
                
        if it >= itmax
           error('Calculo fallido en la curvatura, fin de simulacion!')
        end
    end

normalcen = normalnew;

end
