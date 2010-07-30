% TODO: cuadrar las constantes
% calcula el delta de fuerza debido a bending
function [deltafcurv,deltafbend,sigma] = deltafbending(geom,parms)

excessarea = (geom.s - geom.areaini)/geom.areaini;
options.tolx = 1e-8;
sigma = parms.bending.sigma;
[sigma,valf] = fzero(@zeroalpha,sigma,options);
% Quitar comentario para mostrar resultados de sigma y exceso de area
sigma
excessarea
deltafcurv = (sigma*parms.rkcurv).*(2.*geom.curv);

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

%     figure(6);
%     grafscfld(geom,geom.lapcurv);
%     axis equal; view(90,0); xlabel('x1'); ylabel('x2'); zlabel('x3'); colorbar;
%     getframe; title('laplace curv');
%     figure(7);
%     grafscfld(geom,deltafcurv);
%     axis equal; view(90,0); xlabel('x1'); ylabel('x2'); zlabel('x3'); colorbar;
%     getframe; title('Tension');
%     figure(8);
%     grafscfld(geom,geom.curv);
%     axis equal; view(90,0); xlabel('x1'); ylabel('x2'); zlabel('x3'); colorbar;
%     getframe; title('Curvatura Media');
%     figure(9);
%     grafscfld(geom,geom.Kg);
%     axis equal; view(90,0); xlabel('x1'); ylabel('x2'); zlabel('x3'); colorbar;
%     getframe; title('Curvatura Gaussiana');
%     figure(10)
%     grafscfld(geom,deltafbend);
%     axis equal; view(90,0); xlabel('x1'); ylabel('x2'); zlabel('x3'); colorbar;
%     getframe; title('Bending');

return

function zeroa = zeroalpha(sigmam)
    % Funcion para calcular el esfuerzo usando la expresion
    % Completa dada por Evans and Rawicz, 1990; Rawicz et al., 2000
    zeroa = log(1 + parms.bending.c*sigmam*geom.s)/(8*pi*parms.bending.kbar)...
     + sigmam/parms.bending.kext - excessarea;
end

end
