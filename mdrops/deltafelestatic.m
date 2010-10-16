% Calcula el deltaf debido a la interaccion electroestatica entre la vesicula y
% el vidrio.
function rdeltafelestat = deltafelestatic(geom,parms)

rkelestat = parms.rkelestat;
l = parms.elestat.l;
psi1 = parms.elestat.psi1;
psi2 = parms.elestat.psi2;
h = centroide(geom);
h = h(3);
elest = [0 0 -1];

fuerzaelest = rkelestat*...
   ((2*psi1*psi2*exp(-l*h)+(psi1^2+psi2^2)*exp(-2*l*h))/(1-exp(-2*l*h)))/geom.s;

rdeltafelestat = sum(geom.normal.*repmat(elest,geom.numnodes,1),2).*fuerzaelest;

% rdeltafelestat = -sum(geom.normal.*repmat(ele,geom.numnodes,1),2)...
%    *fuerza/geom.s;