% calcula e deltaf en la interfase de una gota debido a gravedad
function rdeltafgrav = deltafgrav(rnodes,rkgrav)

rnumnodes = size(rnodes,1);

% Vector Unitario de Gravedad
rg=[0 0 1];  
rdeltafgrav = sum(rnodes.*repmat(rg,rnumnodes,1),2).*rkgrav;