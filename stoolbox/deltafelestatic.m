% Calcula el deltaf debido a la interaccion electroestatica entre la vesicula y
% el vidrio.
function [rdeltafelestat,fuerzaelest] = deltafelestatic(geom,parms)

rkelestat = parms.rkelestat;
l = parms.elestat.l;
h = geom.nodes(:,3);
relest = [0 0 -1];

fuerzaelest = l*rkelestat*exp(-l.*h);

% Fuerza volumetrica
% rdeltafelestat = (geom.nodes*relest').*(fuerzaelest/geom.vol);

% Fuerza superficial
rdeltafelestat = (geom.normal*relest').*fuerzaelest;

fuerzaelest = max(fuerzaelest);
