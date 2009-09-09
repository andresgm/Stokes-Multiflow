% devuelve el valor de scfld(t + deltat) mediante el metodo theta:
% theta = 0 Explicit Euler
% theta = 0.5 Crank - Nicholson Semi Implicit
% theta = 1 Full Implicit
% ajimat: matriz global de elementos o volumenes finitos.

function scfld = thetamethod(ajimat,deltat,theta,scfld)

numnodes = size(scfld,1);
kt = (deltat.*(1 - theta));

if theta == 0
    % euler explicit
    scfld = scfld + deltat.*(ajimat*scfld);
else
    % Si theta = 1/2 Crank-Nicolson semi implicito
    % Si theta = 1 full implicito
    dji = sparse(1:1:numnodes,1:1:numnodes,ones(numnodes,1));
    matglob = (dji - (deltat.*theta).*ajimat);
    vector = kt.*(ajimat*scfld) + scfld;
    scfld = matglob\vector;
end
