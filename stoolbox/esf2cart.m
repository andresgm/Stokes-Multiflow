% transforma un arreglo de coordenadas cartesianas en esfericas
% TODO: SUJETO A VERIFICACION
function [x,y,z] = esf2cart(r,teta,phi)

x = r.*cos(teta).*sin(phi);
y = r.*sin(teta).*sin(phi);
z = r.*cos(phi);
