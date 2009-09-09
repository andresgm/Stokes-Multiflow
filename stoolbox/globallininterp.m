function fldinterp = globallininterp(geom,fld,zitas)

[i,jj] = size(geom.elements);
[ii,j] = size(fld);
if nargin < 2
    % por defecto en el centro del triangulo
    zitas = [1/3;1/3];
    nzitas = 1;
else
    nzitas = size(zitas,2);
end
   
fldn1 = fld(geom.elements(:,1),:);
fldn2 = fld(geom.elements(:,2),:);
fldn3 = fld(geom.elements(:,3),:);

fldinterp = zeros(i,j,nzitas);
for k = 1:nzitas
    fldinterp(:,:,k) = lininterp(fldn1,fldn2,fldn3,zitas(:,k));
end
