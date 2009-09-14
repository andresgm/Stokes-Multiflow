% calcula las area de los triangulos definidos por un arreglo de nodos
function triarea = trianglearea(nodo1,nodo2,nodo3)

x1 = nodo1(:,1);
x2 = nodo2(:,1);
x3 = nodo3(:,1);
y1 = nodo1(:,2);
y2 = nodo2(:,2);
y3 = nodo3(:,2);
z1 = nodo1(:,3);
z2 = nodo2(:,3);
z3 = nodo3(:,3);

det1q = ((z2 - z3).*y1 - (y2 - y3).*z1 + y2.*z3 - y3.*z2).^2;

det2q = (z1.*(x2 - x3) - x1.*(z2 - z3) + z2.*x3 - z3.*x2).^2;

det3q = (x1.*(y2 - y3) - y1.*(x2 - x3) + x2.*y3 - x3.*y2).^2;

triarea = 0.5.*(det1q + det2q + det3q).^0.5;
