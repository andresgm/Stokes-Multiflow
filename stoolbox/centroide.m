% calcula el centroide de una malla cerrada
% TODO generalizar para mas "gotas"
function xc = centroide(geom)

if isfield(geom,'element2node') ~= 1
    % calcule la conectividad de elmeentos vecino
    geom.element2node = element2node(geom.elements);
end

if isfield(geom,'normal') == 1
    normnode = geom.normal;
else
    normnode = normal(geom);
end

if isfield(geom,'dsi') == 1
    dsi = geom.dsi;
else
    [s,ds,dsi] = areas(geom); 
end
numnodes = size(geom.nodes,1);

xdotn = sum(geom.nodes.*geom.normal,2);
volcalc = sum(xdotn.*dsi)/3;

    % Calcule el centroide xc Yc Zc
x1 = [geom.nodes(:,1).^2 zeros(numnodes,2)];
x2 = [zeros(numnodes,1) geom.nodes(:,2).^2 zeros(numnodes,1)];
x3 = [zeros(numnodes,2) geom.nodes(:,3).^2];

xc = zeros(1,3);
x1dotn = sum(x1.*normnode,2);
x2dotn = sum(x2.*normnode,2);
x3dotn = sum(x3.*normnode,2);
xc(1) = sum(x1dotn.*dsi);
xc(2) = sum(x2dotn.*dsi);
xc(3) = sum(x3dotn.*dsi);

xc = xc./(2*volcalc);
