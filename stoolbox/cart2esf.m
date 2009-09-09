% transforma un arreglo de coordenadas cartesianas en esfericas
% [r,teta,phi] = cart2esf(nodos)
% r = normesp(nodos)
% teta [0 2pi];
% phi [0 pi];
% TODO: SUJETO A VERIFICACION
function [r,teta,phi] = cart2esf(nodos)

r = normesp(nodos);

teta = atan2(nodos(:,2),nodos(:,1));

q =  nodos(:,2) < 0;

% tercer y cuarto cuadrante
teta(q) = 2.*pi + teta(q);


% % tercer cuadrante
% teta(q3) = 2.*pi + teta(q3);
% 
% % cuarto cuadrante
% teta(q4) =  2.*pi + teta(q4);

phi = acos(nodos(:,3)./r);
