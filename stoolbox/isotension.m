% Funcion para el calculo de la tension isotropica necesaria para la
% conservacion local del area. Esta funcion esta basada en lo propuesto
% en Evans y Skalak 1980 para la componente de tension isotropica.
function [isotens] = isotension(geom,kaes,mues)

% Entrada:
% geom: estructura que contiene entre otras variables las matrices de
%       coordenadas y coordinacion de la malla.

% El siguiente algoritmo sigue la metodologia propuesta en Charrier et al.
% Journal of Strain Analysis Vol 24 No 2 1983

tensionelas = zeros(3,geom.numnodes);

for i=1:size(geom.elements,1)
    
    A = geom.shapeA(:,i);
    B = geom.shapeB(:,i);
    
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
    vE4 = (vXk-vXi)./norm(vXk-vXi);
    vE3 = cross(vE1,vE4)/norm(cross(vE1,vE4));
    vE2 = cross(vE3,vE1);
    
    R = [vE1;vE2;vE3];
    
    vxil = geom.refrot(:,1,i);
    vxjl = geom.refrot(:,2,i);
    vxkl = geom.refrot(:,3,i);
    
    vXil = [0,0,0];
    vXjl = R*(vXj-vXi)';
    vXkl = R*(vXk-vXi)';
    
    u = [vXil(1)-vxil(1),vXjl(1)-vxjl(1),vXkl(1)-vxkl(1)]';
    v = [vXil(2)-vxil(2),vXjl(2)-vxjl(2),vXkl(2)-vxkl(2)]';
    
    g11 = 1 + 2*u'*A + (u'*A)*(A'*u) + (v'*A)*(A'*v);
    g22 = 1 + 2*v'*B + (v'*B)*(B'*v) + (u'*B)*(B'*u);
    g12 = u'*B + (u'*B)*(A'*u) + v'*A + (v'*B)*(A'*v);
    
    l1 = sqrt((g11+g22+sqrt((g11-g22)^2+4*g12^2))/2);
    l2 = sqrt((g11+g22-sqrt((g11-g22)^2+4*g12^2))/2);
    
%     isotens(geom.elements(i,1)) = isotens(geom.elements(i,1))+(l1*l2-1);
%     isotens(geom.elements(i,2)) = isotens(geom.elements(i,2))+(l1*l2-1);
%     isotens(geom.elements(i,3)) = isotens(geom.elements(i,3))+(l1*l2-1);
    dg11dui = 2*A(1)*(1+u'*A);
    dg11duj = 2*A(2)*(1+u'*A);
    dg11duk = 2*A(3)*(1+u'*A);
    dg11dvi = 2*A(1)*(v'*A);
    dg11dvj = 2*A(2)*(v'*A);
    dg11dvk = 2*A(3)*(v'*A);
    
    dg22dui = 2*B(1)*(u'*B);
    dg22duj = 2*B(2)*(u'*B);
    dg22duk = 2*B(3)*(u'*B);
    dg22dvi = 2*B(1)*(1+v'*B);
    dg22dvj = 2*B(2)*(1+v'*B);
    dg22dvk = 2*B(3)*(1+v'*B);
    
    dg12dui = B(1)*(1+u'*A) + A(1)*(u'*B);
    dg12duj = B(2)*(1+u'*A) + A(2)*(u'*B);
    dg12duk = B(3)*(1+u'*A) + A(3)*(u'*B);
    dg12dvi = A(1)*(1+v'*B) + B(1)*(v'*A);
    dg12dvj = A(2)*(1+v'*B) + B(2)*(v'*A);
    dg12dvk = A(3)*(1+v'*B) + B(3)*(v'*A);
    
    cuad = (g11-g22)^2+4*g12^2;
    
    if cuad <= 1e-6
        radical = 0;
    else
        radical = cuad^(-.5);
    end
    
    g11mg22 = g11-g22;
    
    dl1dui = (dg11dui + dg22dui +...
        radical*((g11mg22)*(dg11dui-dg22dui)+4*g12*dg12dui))/(4*l1);
    dl1duj = (dg11duj + dg22duj +...
        radical*((g11mg22)*(dg11duj-dg22duj)+4*g12*dg12duj))/(4*l1);
    dl1duk = (dg11duk + dg22duk +...
        radical*((g11mg22)*(dg11duk-dg22duk)+4*g12*dg12duk))/(4*l1);
    dl1dvi = (dg11dvi + dg22dvi +...
        radical*((g11mg22)*(dg11dvi-dg22dvi)+4*g12*dg12dvi))/(4*l1);
    dl1dvj = (dg11dvj + dg22dvj +...
        radical*((g11mg22)*(dg11dvj-dg22dvj)+4*g12*dg12dvj))/(4*l1);
    dl1dvk = (dg11dvk + dg22dvk +...
        radical*((g11mg22)*(dg11dvk-dg22dvk)+4*g12*dg12dvk))/(4*l1);
    dl2dui = (dg11dui + dg22dui -...
        radical*((g11mg22)*(dg11dui-dg22dui)+4*g12*dg12dui))/(4*l2);
    dl2duj = (dg11duj + dg22duj -...
        radical*((g11mg22)*(dg11duj-dg22duj)+4*g12*dg12duj))/(4*l2);
    dl2duk = (dg11duk + dg22duk -...
        radical*((g11mg22)*(dg11duk-dg22duk)+4*g12*dg12duk))/(4*l2);
    dl2dvi = (dg11dvi + dg22dvi -...
        radical*((g11mg22)*(dg11dvi-dg22dvi)+4*g12*dg12dvi))/(4*l2);
    dl2dvj = (dg11dvj + dg22dvj -...
        radical*((g11mg22)*(dg11dvj-dg22dvj)+4*g12*dg12dvj))/(4*l2);
    dl2dvk = (dg11dvk + dg22dvk -...
        radical*((g11mg22)*(dg11dvk-dg22dvk)+4*g12*dg12dvk))/(4*l2);
    
    % En las lineas a continuacion kaes y mues son los parametros 
    % del modelo de Evans y Skalak.
    % strain energy function from evans and skalak
    % mech and therm of biomembranes 1980 pg 98 (4.8.1)
    dwdl1 = 2*kaes*(l1*l2-1)*l2 - mues*(l2^2)/l1;
    dwdl2 = 2*kaes*(l1*l2-1)*l1 - mues*(l1^2)/l2;
    
    dwdui = dwdl1*dl1dui + dwdl2*dl2dui;
    dwduj = dwdl1*dl1duj + dwdl2*dl2duj;
    dwduk = dwdl1*dl1duk + dwdl2*dl2duk;
    dwdvi = dwdl1*dl1dvi + dwdl2*dl2dvi;
    dwdvj = dwdl1*dl1dvj + dwdl2*dl2dvj;
    dwdvk = dwdl1*dl1dvk + dwdl2*dl2dvk;
    
    fxi = dwdui*geom.dsref(i);
    fyi = dwdvi*geom.dsref(i);
    fzi = 0.0;
    fxj = dwduj*geom.dsref(i);
	fyj = dwdvj*geom.dsref(i);
	fzj = 0.0;
	fxk = dwduk*geom.dsref(i);
	fyk = dwdvk*geom.dsref(i);
	fzk = 0.0;
    
    tensionelas(:,geom.elements(i,1)) = ...
        tensionelas(:,geom.elements(i,1)) - (R'*[fxi;fyi;fzi]);
    tensionelas(:,geom.elements(i,2)) = ...
        tensionelas(:,geom.elements(i,2)) - (R'*[fxj;fyj;fzj]);
    tensionelas(:,geom.elements(i,3)) = ...
        tensionelas(:,geom.elements(i,3)) - (R'*[fxk;fyk;fzk]);
end

isotens = tensionelas;