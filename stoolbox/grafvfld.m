% grafica el campo vectorial vfld definido en los nodos de una malla
% triangular definida por struct

function grafvfld(struct,vfld,options)

quiver3(struct.nodes(:,1),struct.nodes(:,2),struct.nodes(:,3),...
    vfld(:,1),vfld(:,2),vfld(:,3))