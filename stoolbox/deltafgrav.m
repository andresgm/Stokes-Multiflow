% calcula e deltaf en la interfase de una gota debido a gravedad
function rdeltafgrav = deltafgrav(geom,rkgrav)

% Vector Unitario de Gravedad
rg=[0 0 1];
h = centroide(geom);
h = h(3);

fuerzag = h*rkgrav;

rdeltafgrav = sum(geom.normal.*repmat(rg,geom.numnodes,1),2).*fuerzag;

% rdeltafgrav = sum(rnodes.*repmat(rg,geom.numnodes,1),2).*rkgrav;