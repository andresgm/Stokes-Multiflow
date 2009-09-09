% TODO: cuadrar las constantes
% calcula el delta de fuerza debido a bending
function [deltafcurv,deltafbend,sigma] = deltafbending(geom,parms)

excessarea = (geom.s - geom.areaini)/geom.areaini;
options.tolx = 1e-8;
[sigma,valf] = fzero(@zeroalpha,parms.bending.sigma,options);
% Quitar comentario para mostrar resultados de sigma y exceso de area
% Sigma 
% excessarea
parms.bending.sigma = sigma;
deltafcurv = (parms.bending.sigma*parms.rkcurv).*(2.*geom.curv);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % ORIGINAL SIN CURVATURA GAUSSIANA
% deltafbend = -(parms.rkbend).*...
%            (4.*geom.curv.^3 + 2.*geom.lapcurv);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NUEVO CON CURVATURA GAUSSIANA
deltafbend = -(parms.rkbend).*...
           (4.*geom.curv.^3 + 2.*geom.lapcurv - 4.*geom.Kg.*geom.curv);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return

function zeroa = zeroalpha(sigmam)
    % Funcion para calcular el esfuerzo usando la expresion
    % Completa dada por Evans and Rawicz, 1990; Rawicz et al., 2000
    zeroa = log(1 + parms.bending.c*sigmam*geom.s)/(8*pi*parms.bending.kbar)...
     + sigmam*((1.66e-9)*sqrt(parms.bending.kbar))/sqrt(parms.g0)...
     - excessarea;
end

end
