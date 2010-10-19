% calcula e deltaf en la interfase de una gota debido a gravedad
function [rdeltafgrav,fuerzagrav] = deltafgrav(geom,rkgrav)

% Vector Unitario de Gravedad
rg=[0 0 1];

rdeltafgrav = sum(geom.nodes.*repmat(rg,geom.numnodes,1),2).*rkgrav;

fuerzagrav = sum(sum(geom.normal.*repmat(rg,geom.numnodes,1),2)...
    .*rdeltafgrav);