function [abase,abasectr,metricg,metricgctr,jacmet,jacmetctr] = isobase(geom)
% Function to calculate the base vector components in the isoparametric
% triangular element from the global coordinates.

nodes = geom.nodes;
elements = geom.elements;
numelements = geom.numelements;
abase = zeros(3,2,numelements);
abasectr = zeros(3,2,numelements);
metricg = zeros(3,numelements);
metricgctr = zeros(3,numelements);
jacmet = zeros(numelements,1);
jacmetctr = zeros(numelements,1);

for i = 1:numelements
    dN = p1dbasis(0, 0);
    abase(:,1,i) = dN(1,1)*nodes(elements(i,1),:)...
                 + dN(2,1)*nodes(elements(i,2),:)...
                 + dN(3,1)*nodes(elements(i,3),:);
    abase(:,2,i) = dN(1,2)*nodes(elements(i,1),:)...
                 + dN(2,2)*nodes(elements(i,2),:)...
                 + dN(3,2)*nodes(elements(i,3),:);
    % Normal a los vectore a1 y a2
%     a3 = cross(abase(:,1,i),abase(:,2,i));
%     a3 = a3/norm(a3);
    % Normal a la superficie
    a3 = geom.normal(geom.elements(i,1),:);
    ca2a3 = cross(abase(:,2,i),a3);
    v = dot(abase(:,1,i),ca2a3);
    abasectr(:,1,i) = ca2a3/v;
    abasectr(:,2,i) = cross(a3,abase(:,1,i))/v;
    metricg(1,i) = dot(abase(:,1,i),abase(:,1,i));
    metricg(2,i) = dot(abase(:,2,i),abase(:,2,i));
    metricg(3,i) = dot(abase(:,1,i),abase(:,2,i));
    jacmet(i) = metricg(1,i)*metricg(2,i)-metricg(3,i)*metricg(3,i);
    metricgctr(1,i) = dot(abasectr(:,1,i),abasectr(:,1,i));
    metricgctr(2,i) = dot(abasectr(:,2,i),abasectr(:,2,i));
    metricgctr(3,i) = dot(abasectr(:,1,i),abasectr(:,2,i));
    jacmetctr(i) = ...
        metricgctr(1,i)*metricgctr(2,i)-metricgctr(3,i)*metricgctr(3,i);
end