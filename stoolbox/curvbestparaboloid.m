% Rutina Calculo de la Normal y Curvatura media por Best Paraboloid Fitting
% Tomado de Zinchenko 1997

function [curvatnodes,normalatnodes] = curvbestparaboloid(geom,bestparopt)

% recupere las variables
nodes = geom.nodes;
numnodes = size(geom.nodes,1);

if nargin < 2
    % no hay entradas de opciones - opciones por defecto
    % tolerancia
    tol = 1e-8;
    % maximo numero de iteraciones
    itmax = 10;
else
    tol = bestparopt.tol;
    itmax = bestparopt.itmax;
end

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

% Defina espacio en memoria
curvatnodes = zeros(numnodes,1);
normalatnodes = zeros(numnodes,3);

% ORIGINAL
for i=1:numnodes
    % extraiga los indices de los nodos vertices conectados al i
    nodesadj = nodecon2node{i};
    numnodesadj = size(nodesadj,1);
    % extraiga las coordenadas de los nodos conectados al is
    xglobnodes = nodes(nodesadj,:);
    if numnodesadj == 4
        % Genere un par de nodos adicionales
        xglobnodes(5,:) = (xglobnodes(2,:) - xglobnodes(1,:)).*0.5 + xglobnodes(1,:);
        xglobnodes(6,:) = (xglobnodes(3,:) - xglobnodes(2,:)).*0.5 + xglobnodes(2,:);
        numnodesadj = 6;
    end
    
    normalnew = norm_i(i,:);
    normalant = zeros(1,3);
    
    xglobnodes = xglobnodes - repmat(nodes(i,:),numnodesadj,1); % Origen coordenadas globales
																	% en nodo i    
    it = 0;                                                                    
    while normesp(normalnew - normalant) > tol
        it = it + 1;
    % xlocnodes = localcoord(nodesadj);
	% Defina la base local en el i a partir de la coordenada del nodo 1
    vectp = nodes(nodesadj(1),:) - nodes(i,:); % esto ya se hizo y solo hay que hacerlo una vez
    temp1 = vectp - sum(normalnew.*vectp,2).*normalnew;
    t1 = temp1./normesp(temp1);
    t2 = cross(normalnew,t1);
    % Determine las coordenadas locales de los nodos adyacentes a la base
    % ubicada en el i
    xlocnodes = [sum(xglobnodes.*repmat(t1,numnodesadj,1),2)...
                 sum(xglobnodes.*repmat(t2,numnodesadj,1),2)...
                 sum(xglobnodes.*repmat(normalnew,numnodesadj,1),2)];
    % invoque armaje y solucion del sistema lineal asociado al i
    alfa = linearassembly(xlocnodes);
    
    % calculo de la normal nueva en la base local
        normalnewl = [-alfa(1) -alfa(2) 1]./(1 + alfa(1)^2 + alfa(2)^2)^0.5; 
    % transforme la normal a la base Global
        normalant = normalnew;
        normalnew = ([t1' t2' normalnew']*normalnewl')';
    % calculo de la curvatura media
        curvnew = -(alfa(3) + alfa(5));
        
        if it >= itmax
           disp('no se alcanzo la convergencia, calculo fallido') 
           error('calculo fallido en la curvatura fin de simulacion')
        end
    end
        curvatnodes(i) = curvnew;
        normalatnodes(i,:) = normalnew./normesp(normalnew);
end


function alfa = linearassembly(xlocnodes)
    
% calculo de ri = (xi^2 + yi^2 + zi^2)

ri = sum(xlocnodes.^2,2);

% coeficientes de la matriz
xi = xlocnodes(:,1);
xiq = xlocnodes(:,1).^2;
yi = xlocnodes(:,2);
yiq = xlocnodes(:,2).^2;
zi = xlocnodes(:,3);

% TODO: Por que no se usa??? OJO REVIZAR 
% ziq = xlocnodes(:,3).^2;

matcoef = [sum(xiq./ri,1) sum(yi.*xi./ri,1) sum(xiq.*xi./ri,1) sum(xiq.*yi./ri,1) sum(yiq.*xi./ri,1);...
           sum(xi.*yi./ri,1) sum(yiq./ri,1) sum(xiq.*yi./ri,1) sum(xi.*yiq./ri,1) sum(yiq.*yi./ri,1);...
           sum(xiq.*xi./ri,1) sum(yi.*xiq./ri,1) sum(xiq.*xiq./ri,1) sum(xiq.*xi.*yi./ri,1) sum(yiq.*xiq./ri,1);...
           sum(xiq.*yi./ri,1) sum(xi.*yiq./ri,1) sum(xiq.*xi.*yi./ri,1) sum(xiq.*yiq./ri,1) sum(yiq.*yi.*xi./ri,1);...
           sum(yiq.*xi./ri,1) sum(yi.*yiq./ri,1) sum(xiq.*yiq./ri,1) sum(xi.*yiq.*yi./ri,1) sum(yiq.*yiq./ri,1)];
vectcoef = [sum(zi.*xi./ri,1); sum(zi.*yi./ri,1); sum(zi.*xiq./ri,1); sum(zi.*xi.*yi./ri,1); sum(zi.*yiq./ri,1)];

alfa = matcoef\vectcoef;
           
% % Opcion de modificar
% for i=1:numnodes
%     % extraiga los indices de los nodos vertices conectados al i
%     nodesadj = nodecon2node{i};
%     % extraiga las coordenadas de los nodos conectados al is
%     xglobnodes = nodes(nodesadj,:);
%     if size(nodesadj,1) == 4
%         % Genere un par de nodos adicionales
%         xglobnodes(5,:) = (xglobnodes(2,:) - xglobnodes(1,:)).*0.5 + xglobnodes(1,:);
%         xglobnodes(6,:) = (xglobnodes(3,:) - xglobnodes(2,:)).*0.5 + xglobnodes(2,:);
%         numnodesadj = 6;
%     end
%        
%     % calcule las posiciones locales
%     % extraiga las coordenadas de los nodos conectados al is
%     xglobnodes = geom.nodes(nodesadj,:)';
%     nodecen = geom.nodes(i,:)';
% 
%     deltax = zeros(3,size(nodesadj,1));
%     for k=1:size(nodesadj,1)
%         temp1 = xglobnodes(:,k) - nodecen;
%         deltax(:,k) = temp1;
%     end
%     
%     normalnew = norm_i(i,:);
%     normalant = zeros(1,3);
%     
%     it = 0;                                                                    
%     while normesp(normalnew - normalant) > tol
%         it = it + 1;
%         % matriz de rotaciones
%         r3 = normalnew;
%         pga = projtensor(normalnew);
%         r1 = pga(:,1)';
%         r1 = univect(r1);
%         r2 = crossv(r3,r1);
%             
%         mr = [r1;r2;r3];
%         
%         xlocnodes = zeros(3,size(deltax,2));
%         for k=1:size(deltax,2)
%             xlocnodes(:,k) = mr*deltax(:,k);
%         end
%         xlocnodes = xlocnodes';
%         
%     % invoque armaje y solucion del sistema lineal asociado al i
%         alfa = linearassembly(xlocnodes);
%     
%     % calculo de la normal nueva en la base local
%         normalnewl = [-alfa(1) -alfa(2) 1]'./(1 + alfa(1)^2 + alfa(2)^2)^0.5; 
%         normalant = normalnew;
%         % transforme la normal a la base Global
%         normalnew = (mr'*normalnewl)';
%         
%     % calculo de la curvatura media
%         curvnew = -(alfa(3) + alfa(5));
%         
%         if it >= itmax
%            disp('no se alcanzo la convergencia, calculo fallido') 
%            error('calculo fallido en la curvatura fin de simulacion')
%         end
%     end
%         curvatnodes(i) = curvnew;
%         normalatnodes(i,:) = univect(normalnew);
% end


