% calcula el volumende un cuerpo encerrado por una malla superficial
function volumen = volume(geom,parfield)

if nargin == 2
    if isfield(parfield,'jacomp') ~= 1
        % no hay metrica de transformacion calculela
        jacomp = metrictrans(geom,[1/3;1/3]);
    else
        jacomp = parfield.jacomp;
    end

    if isfield(parfield,'normal') ~= 1
        % no hay vector normal calculelo
        normalnode = normal(geom,jacomp);
    else
        normalnode = parfield.normal;
    end

    if isfield(parfield,'dsi') ~= 1
        % no hay areas superficiales definidas calculelas
        [s,ds,dsi] = areas(geom);
    else
        dsi = parfield.dsi;
    end
    
elseif nargin == 1
    % no hay metrica de transformacion calculela
        jacomp = metrictrans(geom,[1/3;1/3]);
        
    if isfield(geom,'normal') ~= 1
        % no hay vector normal calculelo
        normalnode = normal(geom,jacomp);
    else
        normalnode = geom.normal;
    end

    if isfield(geom,'dsi') ~= 1
        % no hay areas superficiales definidas calculelas
        [s,ds,dsi] = areas(geom);
    else
        dsi = geom.dsi;
    end
end

% calculo del producto punto <xnode,normalatnodes>
xdotn = sum(geom.nodes.*normalnode,2);
volumen = sum(xdotn.*dsi)/3;
end


