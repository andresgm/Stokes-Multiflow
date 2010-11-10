% TODO: cuadrar las constantes de la fuerza del campo electrico
function deltafelec = deltafelectric(geom,parms)

% Vector Unitario del campo electrico
ele = [0 0 1];

% Calculo de E0*<ele,n>
deltafelec = (geom.normal*ele')*(parms.rkelect);

