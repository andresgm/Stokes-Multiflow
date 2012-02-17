% calcula el integrando del particle stress tensor sobre la superficie 
% definida por nodes ec 86 PhD Tesis Kennedy
% deltaf: array de dim(numnodes,3) del delta de fuerza en la interfase
% nodes: array de dim(numnodes,3) de los nodos de la interfase
% velnode: array de dim(numnodes,3) de la velocidad en la interfase
% normalv: array de dim(numnodes,3) de la normal en la interfase
% const: const(1): viscosidad dinamica del flujo
% const(2): lamda


function strtensornode = strtensorreologia(deltaf,nodes,velnode,normalv,const)

kappa = (1-const)/(1+const);
strtensornode = [];



for i = 1:3
    for j = 1:3
        strtensornode(i,j,:) = nodes(:,j).*deltaf(:,i) - kappa.*...
            (velnode(:,i).*normalv(:,j) + velnode(:,j).*normalv(:,i));
    end
end

% for i = 1:3
%     for j = 1:3
%         strtensornode(i,j,:) = nodes(:,j).*deltaf(:,i) ...
%             - kappa.*...
%             (velnode(:,i).*normalv(:,j) + velnode(:,j).*normalv(:,i));
%     end
% end