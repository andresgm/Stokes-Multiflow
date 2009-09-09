% calcula el integrando del particle stress tensor sobre la superficie 
% definida por nodes ec 86 PhD Tesis Kennedy
% deltaf: array de dim(numnodes,3) del delta de fuerza en la interfase
% nodes: array de dim(numnodes,3) de los nodos de la interfase
% velnode: array de dim(numnodes,3) de la velocidad en la interfase
% normalv: array de dim(numnodes,3) de la normal en la interfase
% const: const(1): viscosidad dinamica del flujo
% const(2): lamda
% TODO: adimensionalizar el integrando (esta en terminos de la viscosidad)

function strtensornode = partstrtensor(deltaf,nodes,velnode,normalv,const)

miu = const(1);
lamda = const(2);
numnodes = size(nodes,1);
strtensornode = zeros(3,3,numnodes);


for i = 1:3
    for j = 1:3
        strtensornode(i,j,:) = nodes(:,j).*deltaf(:,i) ...
            - miu.*(1 - lamda).*...
            (velnode(:,i).*normalv(:,j) + velnode(:,j).*normalv(:,i));
    end
end

