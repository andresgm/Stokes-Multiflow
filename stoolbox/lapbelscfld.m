% calcula el laplace beltrami de un campo escalar scfld
function laplacebelscfld = lapbelscfld(geom,scfld,l)

if nargin < 3
    % calcule la matriz del laplace beltrami
    l = laplacebeltramimat(geom);
end

laplacebelscfld = l*scfld;
