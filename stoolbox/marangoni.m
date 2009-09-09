% Calcula los esfuerzos de Marangoni debido a la concentracion local
% surfactantes. Basado en Bazhlekov "Numerical investigation of the effect 
% of insoluble surfactants on drop deformation and breakup in simple shear 
% ?ow".
% El calculo es realizado sobre una malla triangular struct
% gamma: concentracion de sulfactante en cada nodo de la malla.
% struct: structura de la malla.
% sigmaparms: parametros de la ecuacion de estado de la tension superficial
% (sigma) (Ecuacion 8)

function [maranstress,sigmavar] = marangoni(struct,gamma,sigmaparms)

% Dada la concentracion gamma en la interfase calcule sigma
sigmavar = sigmamodel(gamma,sigmaparms);

% Calcule el gradiente superficial grad_s(sigma) EN LOS NODOS
grad_ssigma = grad_s(struct,sigmavar);

% Esfuerzos de marangoni
gradsnorm = repmat(sum(grad_ssigma.*struct.normal,2),[1 3]).*struct.normal;
grad_ssigma = grad_ssigma - gradsnorm;

maranstress = -grad_ssigma;
%=======
%% maranstress = -grad_ssigma;
%gradsigmagamma = sigmaparms.e.*sigmaparms.x.*log(1 - sigmaparms.x.*gamma)./((1 - sigmaparms.x.*gamma).*(1 + log(1  - sigmaparms.x)));
%grad_sgamma = grad_s(struct,gamma);
%maranstress = -repmat(gradsigmagamma,[1 3]).*grad_sgamma;
%>>>>>>> .r75

% % Esfuerzos de Marangoni (LapBel(sigma))
% maranstress = - laplacebeltrami(struct,sigmavar);


function sigmavar = sigmamodel(gammav,parms)
% Ecuacion (8)
if strcmp(parms.maranmodel,'linear') == 1
    sigmavar = (1 - parms.beta.*gammav)./(1 - parms.beta);
elseif strcmp(parms.maranmodel,'log') == 1
    sigmavar = (1 + parms.e.*log(1 - parms.x.*gammav))./(1 + parms.e.*log(1 - parms.x));
end

% % Ecuacion (8)
% if strcmp(parms.maranmodel,'linear') == 1
%     sigmavar = (1 - parms.beta.*gammav)./(1 - parms.beta);
% elseif strcmp(parms.maranmodel,'log') == 1
%     sigmavar = (1 + parms.e.*log10(1 - parms.x.*gammav))./(1 + parms.e.*log10(1 - parms.x));
% end


