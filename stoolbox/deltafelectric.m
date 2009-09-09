% TODO: cuadrar las constantes de la fuerza del campo electrico
function deltafelec = deltafelectric(geom,parms)

% Vector Unitario del campo electrico
ele = [0 1 0];
% Calculo de E0*<ele,n>
deltafelec = -sum(geom.normal.*repmat(ele,geom.numnodes,1),2)...
            *(parms.rkelect); 
