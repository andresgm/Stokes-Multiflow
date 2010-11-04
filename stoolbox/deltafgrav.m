% calcula e deltaf en la interfase de una gota debido a gravedad
function [rdeltafgrav,fuerzagrav] = deltafgrav(geom,rkgrav)

% Vector Unitario de Gravedad
rg=[0 0 1];

rdeltafgrav = (geom.nodes*rg').*rkgrav;

fuerzagrav = rkgrav*geom.vol;