% interpolacion lineal isoparametrica de un campo escalar o un campo
% vectorial en elementos triangulares lineales
% xi = campo escalar o vectorial en el nodo i local
% Zita = [Zita1 Zita 2]; Zita3 = 1 - Zita1 - Zita2
function x = lininterp(x1,x2,x3,Zita)

x = x1.*Zita(1) + x2.*Zita(2) + x3.*(1 - Zita(1) - Zita(2));
