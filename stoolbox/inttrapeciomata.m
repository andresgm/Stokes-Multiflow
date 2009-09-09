% Integracion mediante el metodo del trapecio sobre una malla superficial
% triangular de 3 nodos. Basado en Rallison caso de entrada arrays
% rintegrand: arreglo tensorial de la forma n x m x numnodes
% rsi: area baricentrica alrededor de cada nodo
% salida integral de dimension 1 x i componentes.

% rintegrand: arreglo tensorial de la forma n x m x numnodes
function rintegral = inttrapeciomata(rdsi,rintegrand)

[n,m,numnodes] = size(rintegrand);

rintegral = sum(rintegrand.*repmat(permute(rdsi,[3 2 1]),[n m 1]),3);