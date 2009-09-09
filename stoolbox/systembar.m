% devuelve la barra '/' '\' segun sea el sistema operativo
function sbar = systembar()

if ispc == 1
    sbar = '\';
elseif isunix == 1
    sbar = '/';
end