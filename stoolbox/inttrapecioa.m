% Integracion mediante el metodo del trapecio sobre una malla superficial
% triangular de 3 nodos. Basado en Rallison caso de entrada arrays
% integrand: escalar o vectorial de dimension m nodos x i componentes
% rsi: area baricentrica alrededor de cada nodo
% salida integral de dimension 1 x i componentes.
function rintegral = inttrapecioa(rdsi,rintegrand)

[numnodes,comp,s] = size(rintegrand);

rintegral = sum(rintegrand.*repmat(rdsi,[1 comp]),1);