% calcula las direcciones principales de deformacion, el tensor inercial  y
% la deformacion de una malla cerrada.

function [inerttensor,def,v,alpha] = dirprindef(geom,opalpha)

numnodes = size(geom.nodes,1);

% calcule el tensor de inercia centroidal
[inerttensor,nodesloc] = itensor(geom);
% calcule las direcciones principales de deformacion y los autovalores
[v,eigval] = eig(inerttensor);
v = v./repmat((normesp(v'))',[3 1]);

nodesb = zeros(geom.numnodes,3);
for i=1:geom.numnodes
   nodesb(i,:) = (v*geom.nodes(i,:)')';
end

nu(1) = max(abs(sum(nodesloc.*repmat(v(:,1)',numnodes,1),2)));
nu(2) = max(abs(sum(nodesloc.*repmat(v(:,2)',numnodes,1),2)));
nu(3) = max(abs(sum(nodesloc.*repmat(v(:,3)',numnodes,1),2)));

% disp('defeformacion (L - B)/(L + B)');
def = (max(nu) - min(nu))/(max(nu) + min(nu));

if nargin > 1
    % calculo del agulo theta de deformacion en el plano x2-x3
    maxdef = v(:,find(nu == max(nu)));
    alpha = atand(maxdef(3)/maxdef(2));
else
    alpha = 'not calc';
end