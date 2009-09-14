% calcula la curvatura media mediante laplace beltrami
% TODO dar mas salidas 
function curv = curvlb(geom,l)

if isfield(geom,'normal') ~= 1
    % no existe calculo de la normal ejecutelo
    normal_i = normal(geom);
else
   normal_i = geom.normal; 
end

if nargin < 2
    % calcule la matriz del laplace beltrami
    l = laplacebeltramimat(geom);
end

kx1 = l*geom.nodes(:,1);
kx2 = l*geom.nodes(:,2);
kx3 = l*geom.nodes(:,3);

% Curvature normal
k = [kx1 kx2 kx3];
% Mean Curvature
curv = normesp(k)./2;

% verifique que la curvatura es positiva
curvnegind = sum(normal_i.*k,2) >= 0;
curv(curvnegind) = -curv(curvnegind);
