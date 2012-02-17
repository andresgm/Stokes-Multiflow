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

        tol = 1e-8;
        itmax = 100;
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

    deltax = xglobnodes-repmat(nodecen,1,size(xglobnodes,2));
%     for k=1:size(nodesadj,1)
%         temp1 = xglobnodes(:,k) - nodecen;
%         deltax(:,k) = temp1;
%     end

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
        % [cmean(i),normalbp(i,:),Kg(i)] = extendedquatric(deltax,norm_i(i,:),tol,itmax);
        [normalbp(i,:),cmean(i),Kg(i)] = bestparaboloid(deltax,norm_i(i,:),tol,itmax);
    end

end

if strcmp(tipo,'single') == 1
    % Si se usa metodo 'single' la normal no se actualiza
    normal_f = norm_i;
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

% function [cmeancen,normalcen,Kg] = extendedquatric(deltax,normalnew,tol,itmax)
% 
%     normalant = [0,0,0];
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
%     % ensamble la matriz del sistema lineal
%         col1 = xlocnodes(:,1).^2;
%         col2 = xlocnodes(:,1).*xlocnodes(:,2);
%         col3 = xlocnodes(:,2).^2;
%         col4 = xlocnodes(:,1);
%         col5 = xlocnodes(:,2);
%         
%         matriz = [col1 col2 col3 col4 col5];
%         zeta = xlocnodes(:,3);
% 
%         % Solucion en minimos cuadrados
%         alfa = matriz\zeta;
% 
%         % extraccion de propiedades
%         % normal en la base local
%         normalant = normalnew;
%         normalnewl = (1/(alfa(4)^2 + alfa(5)^2 + 1)^0.5).*[-alfa(4),-alfa(5),1]';
%         % transforme la normal a la base Global
%         normalnew = (mr'*normalnewl)';
%         % curvatura media
%         cmeancen = -(alfa(1) + alfa(3) + alfa(1)*alfa(5)^2 + alfa(3)*alfa(4)^2 ...
%             - alfa(2)*alfa(4)*alfa(5))/(1 + alfa(4)^2 + alfa(5)^2)^1.5;
%         
%         
%         % curvatura gaussiana
%         Kg= ((4.*alfa(1).*alfa(3) - alfa(2)^2)./(1 + alfa(4)^2 + alfa(5)^2)^2);
%                 
%         if it >= itmax
%            error('Calculo fallido en la curvatura, fin de simulacion!')
%         end
%     end
% normalcen = normalnew;
% 
% end

function [normal,cmean,Kg] = bestparaboloid(nbcoord,normal,tol,itmax)
% Calcula la normal y curvaturas media y gaussiana usando el metodo de ajuste al
% mejor paraboloide usando minimos cuadrados.
    normalant = [0,0,0];
    normalnew = normal;
    it = 0;
    while normesp(normalnew - normalant) > tol
    it = it + 1;
    % matriz de rotaciones
       r3 = normalnew;
       pga = eye(3)-normalnew'*normalnew;
       %r1 = pga(:,1)';
       r1 = (pga*nbcoord(:,1))';
       r1 = univect(r1);
       r2 = crossv(r3,r1);
            
       mr = [r1;r2;r3];
        
       xlocnodes = zeros(3,size(nbcoord,2));
       for k=1:size(nbcoord,2)
           xlocnodes(:,k) = mr*nbcoord(:,k);
       end
       xlocnodes = xlocnodes';
       mata = zeros(5,5);
       vectb = zeros(5,1);
       for i = 1:size(xlocnodes,1)
          xi = xlocnodes(i,1);
          yi = xlocnodes(i,2);
          zi = xlocnodes(i,3);
          dist = xi^2+yi^2+zi^2;
          mata(1,1) = mata(1,1) + xi^4/dist;
          mata(1,2) = mata(1,2) + xi^3*yi/dist;
          mata(1,3) = mata(1,3) + yi^2*xi^2/dist;
          mata(1,4) = mata(1,4) + xi^3/dist;
          mata(1,5) = mata(1,5) + xi^2*yi/dist;
          mata(2,1) = mata(2,1) + xi^3*yi/dist;
          mata(2,2) = mata(2,2) + xi^2*yi^2/dist;
          mata(2,3) = mata(2,3) + xi*yi^3/dist;
          mata(2,4) = mata(2,4) + xi^2*yi/dist;
          mata(2,5) = mata(2,5) + xi*yi^2/dist;
          mata(3,1) = mata(3,1) + xi^2*yi^2/dist;
          mata(3,2) = mata(3,2) + xi*yi^3/dist;
          mata(3,3) = mata(3,3) + yi^4/dist;
          mata(3,4) = mata(3,4) + xi*yi^2/dist;
          mata(3,5) = mata(3,5) + yi^3/dist;
          mata(4,1) = mata(4,1) + xi^3/dist;
          mata(4,2) = mata(4,2) + xi^2*yi/dist;
          mata(4,3) = mata(4,3) + xi*yi^2/dist;
          mata(4,4) = mata(4,4) + xi^2/dist;
          mata(4,5) = mata(4,5) + xi*yi/dist;
          mata(5,1) = mata(5,1) + xi^2*yi/dist;
          mata(5,2) = mata(5,2) + xi*yi^2/dist;
          mata(5,3) = mata(5,3) + yi^3/dist;
          mata(5,4) = mata(5,4) + xi*yi/dist;
          mata(5,5) = mata(5,5) + yi^2/dist;
          vectb(1) = vectb(1) + xi^2*zi/dist;
          vectb(2) = vectb(2) + xi*yi*zi/dist;
          vectb(3) = vectb(3) + yi^2*zi/dist;
          vectb(4) = vectb(4) + xi*zi/dist;
          vectb(5) = vectb(5) + yi*zi/dist;
       end
       paramvec = mata\vectb;
       a = paramvec(1);
       b = paramvec(2);
       c = paramvec(3);
       d = paramvec(4);
       e = paramvec(5);
       normalant = normalnew;
       normalnewl = [-d,-e,1]'./(sqrt(d^2 + e^2 + 1));
       % transforme la normal a la base Global
       normalnew = (mr'*normalnewl)';
       % curvatura media
       cmean = -(a+c+a*e^2+c*d^2-b*d*e)/(1 + d^2 + e^2)^1.5;
       % curvatura gaussiana
       Kg= (4*a*c-b^2)/(1 + d^2 + e^2)^2;
       if it >= itmax
          error('Calculo fallido en la curvatura, fin de simulacion!')
       end
   end
   normal = normalnew;


end