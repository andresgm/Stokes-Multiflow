% grafica el campo escalar scfld definido en los nodos de una malla
% triangular definida por struct

function grafscfld(struct,scfld,options)

trisurf(struct.elements,struct.nodes(:,1),struct.nodes(:,2),...
    struct.nodes(:,3),scfld)