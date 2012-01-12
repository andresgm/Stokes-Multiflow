% Funcion para el calculo de la tension isotropica necesaria para la
% conservacion local del area. Esta funcion esta basada en lo propuesto
% en Evans y Skalak 1980 para la componente de tension isotropica.
function [isotens] = isotension(geom,kaes)

% Entrada:
% geom: estructura que contiene entre otras variables las matrices de
%       coordenadas y coordinacion de la malla.

% El siguiente algoritmo sigue la metodologia propuesta en Charrier et al.
% Journal of Strain Analysis Vol 24 No 2 1983

isotens = zeros(geom.numnodes,1);

for i=1:size(geom.elements,1)
    
    A = geom.shapeA(i,:)';
    B = geom.shapeB(i,:)';
    
    Xi = geom.nodes(geom.elements(i,1),1);
    Yi = geom.nodes(geom.elements(i,1),2);
    Zi = geom.nodes(geom.elements(i,1),3);
    
    vXi = [Xi,Yi,Zi];
    
    Xj = geom.nodes(geom.elements(i,2),1);
    Yj = geom.nodes(geom.elements(i,2),2);
    Zj = geom.nodes(geom.elements(i,2),3);
    
    vXj = [Xj,Yj,Zj];
    
    Xk = geom.nodes(geom.elements(i,3),1);
    Yk = geom.nodes(geom.elements(i,3),2);
    Zk = geom.nodes(geom.elements(i,3),3);
    
    vXk = [Xk,Yk,Zk];
    
    vE1 = (vXj-vXi)./norm(vXj-vXi);
    vEtemp = (vXk-vXi)./norm(vXk-vXi);
    vE3 = cross(vE1,vEtemp);
    vE2 = cross(vE3,vE1);
    
    R = [vE1;vE2;vE3];
    
    vxil = geom.refrot(i,1,:);
    vxjl = geom.refrot(i,2,:);
    vxkl = geom.refrot(i,3,:);
    
    vXil = [0,0,0];
    vXjl = R*(vXj-vXi)';
    vXkl = R*(vXk-vXi)';
    
    u = [vXil(1)-vxil(1),vXjl(1)-vxjl(1),vXkl(1)-vxkl(1)]';
    v = [vXil(2)-vxil(2),vXjl(2)-vxjl(2),vXkl(2)-vxkl(2)]';
    
    g11 = 1 + 2*u'*A + (u'*A)*(A'*u) + (v'*A)*(A'*v);
    g22 = 1 * 2*v'*B + (v'*B)*(B'*v) + (u'*B)*(B'*u);
    g12 = u'*B + (u'*B)*(A'*u) + v'*A + (v'*B)*(A'*v);
    
    l1 = sqrt((g11+g22+sqrt((g11-g22)^2+4*g12^2))/2);
    l2 = sqrt((g11+g22-sqrt((g11-g22)^2+4*g12^2))/2);
    
    isotens(i) = (l1*l2-1);    
end

isotens = isotens*kaes;