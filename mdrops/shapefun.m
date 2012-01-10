% Funcion que calcula los coeficientes de las funciones de forma para elementos
% finitos lineales.
function [shapeA, shapeB, refrot] = shapefun(geom)

refrot = zeros(size(geom.elements,1),3,3);
shapeA = zeros(size(geom.elements,1),3);
shapea = zeros(size(geom.elements,1),3);

for i=1:size(geom.elements,1)
    xi = geom.ref(geom.elements(i,1),1);
    yi = geom.ref(geom.elements(i,1),2);
    zi = geom.ref(geom.elements(i,1),3);
    
    vxi = [xi,yi,zi];
    
    xj = geom.ref(geom.elements(i,2),1);
    yj = geom.ref(geom.elements(i,2),2);
    zj = geom.ref(geom.elements(i,2),3);
    
    vxj = [xj,yj,zj];
    
    xk = geom.ref(geom.elements(i,3),1);
    yk = geom.ref(geom.elements(i,3),2);
    zk = geom.ref(geom.elements(i,3),3);
    
    vxk = [xk,yk,zk];
    
    ve1 = (vxj-vxi)./norm(vxj-vxi);
    vetemp = (vxk-vxi)./norm(vxk-vxi);
    ve3 = cross(ve1,vetemp);
    ve2 = cross(ve3,ve1);
    
    r = [ve1;ve2;ve3];
    
    vxil = [0,0,0];
    vxjl = r*(vxj-vxi)';
    vxkl = r*(vxk-vxi)';
    
    refrot(i,:,:) = [vxil;vxjl';vxkl'];
    
    ai = vxjl(2) - vxkl(2);
    bi = vxkl(1) - vxjl(1);
    ci = vxjl(1)*vxkl(2) - vxkl(1)*vxjl(2);
    Li = ci; % Notese que tango xi como yi son zero en el marco de referencia 
             % rotado.
    
    aj = vxkl(2) - vxil(2);
    bj = vxil(1) - vxkl(1);
    cj = vxkl(1)*vxil(2) - vxil(1)*vxkl(2);
    Lj = aj*vxjl(1)+bj*vxjl(2)+cj;
    
    ak = vxil(2) - vxjl(2);
    bk = vxjl(1) - vxil(1);
    ck = vxil(1)*vxjl(2) - vxjl(1)*vxil(2);
    Lk = ak*vxkl(1)+bk*vxkl(2)+ck;
    
    shapeA(i,:) = [ai/Li,aj/Lj,ak/Lk];
    shapeB(i,:) = [bi/Li,bj/Lj,bk/Lk];
end