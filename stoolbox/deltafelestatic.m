% Calcula el deltaf debido a la interaccion electroestatica entre la vesicula y
% el vidrio.
function [rdeltafelestat,fuerzaelest] = deltafelestatic(geom,parms)

rkelestat = parms.rkelestat;
l = parms.elestat.l;
psi1 = parms.elestat.psi1;
psi2 = parms.elestat.psi2;
h = min(geom.nodes(:,3));
elest = [0 0 -1];

propnorm = geom.normal*elest';
areanorm = sum(abs(propnorm).*geom.dsi);

fuerzaelest = rkelestat*...
   ((2*psi1*psi2*exp(-l*h)+(psi1^2+psi2^2)*exp(-2*l*h))/(1-exp(-2*l*h)));

rdeltafelestat = propnorm*fuerzaelest/areanorm;
