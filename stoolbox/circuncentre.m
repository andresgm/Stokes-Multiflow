% calcular el circuncentro de los triangulos definidos por vertices
% p1,p2,p3
function x = circuncentre(p1,p2,p3)

p1_p2 = p1 - p2;
p2_p3 = p2 - p3;
p1_p3 = p1 - p3;

c1 = (normesp(p2_p3).^2 .* sum(p1_p2.*p1_p3,2));
c2 = 2.*normesp(cross(p1_p2,p2_p3)).^2;
alpha = c1./c2;
c1 = (normesp(p1_p3).^2 .* sum(-p1_p2.*p2_p3,2));
beta = c1./c2;
c1 = (normesp(p1_p2).^2 .* sum(-p1_p3.*(-p2_p3),2));
gamma = c1./c2;

x = repmat(alpha,[1 3]).*p1 + repmat(beta,[1 3]).*p2 + repmat(gamma,[1 3]).*p3;

return
