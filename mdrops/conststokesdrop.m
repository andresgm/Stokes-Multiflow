function parms = conststokesdrop(adim,flow,ca,lamda,g0,e0,ka,kb,kd)

if adim == 1
    % adimensionalizacion de andres gonzalez
    parms.flow = flow;
    parms.w = 0;
    % adimensionalizacion del single layer
    parms.rkextf = 2*ca;
    parms.rksl = 2;
    parms.rkdl = 2*(lamda - 1)/(lamda + 1);

    parms.g0 = g0;   
    parms.lamda = lamda;
    parms.ca = ca;
    % parametros de simulacion
    if ka == 1 
        % curvatura constante
        parms.rkcurv = 1/g0;
    else
        parms.rkcurv = 0;
    end

    if kb == 1
        % gravedad
        parms.rkgrav = 1;
    else
        parms.rkgrav = 0;
    end    
      
    if kd == 1
        % campo electrico
        parms.rkelect = e0/g0;
    else
        parms.rkelect = 0;
    end    
     warning...
        ('Se considera la adimensionalizacion: DeltaF = (Sigma/g0)(2H) + Z - (1/g0)(4H^3 + 2Laps(H))');
elseif adim == 2
    % adimensionalizacion de kennedy - No gravity
    parms.flow = flow;
    parms.w = 0;
    % parametros
    % capilar
    parms.ca = ca;
    parms.g0 = 1;
    % lamda
    parms.lamda = lamda;
    
    % constantes de las integrales
    % constante flujo externo
    parms.rkextf = 2*ca;
    % constante single layer
    parms.rksl = 2;
    % constante double layer
    parms.rkdl = 2*(lamda - 1)/(lamda + 1);
    
    % constantes de delta de fuerza
    % constante curvatura
    if ka == 1
        parms.rkcurv = 1;
    else
        parms.rkcurv = 0;
    end
    
    % constante gravedad
    parms.rkgrav = 0;
    
    if kd == 1
        % campo electrico
        parms.rkelect = e0;
    else
        parms.rkelect = 0;
    end    
    
    warning('Kennedy no usa el termino de gravedad rkgrav = 0 por defecto');
elseif adim == 3
    % adimensionalizacion de bazhlekov - shear
    parms.flow = flow;
    parms.w = 0;
    % parametros
    % capilar
    parms.ca = ca;
    parms.g0 = g0;  
    % lamda
    parms.lamda = lamda;
    
    % constantes de las integrales
    % constante flujo externo
    parms.rkextf = 2/(lamda + 1);
    % constante single layer
    parms.rksl = 2/(lamda + 1);
    % constante double layer
    parms.rkdl = 2*(lamda - 1)/(lamda + 1);
    
    % constantes de delta de fuerza
    % constante curvatura
    if ka == 1
        parms.rkcurv = 1/ca;
    else
        parms.rkcurv = 0;
    end
    
    % constante gravedad
    parms.rkgrav = 0;
    
    if kd == 1
        % campo electrico
        parms.rkelect = e0;
    else
        parms.rkelect = 0;
    end    
    
    warning('Adimensionalizacion de Bazhlekov Shear no usa el termino de gravedad rkgrav = 0 por defecto');    
elseif adim == 4
    % adimensionalizacion de bazhlekov - gravity
    parms.flow = flow;
    parms.w = 0;
    % parametros
    % capilar
    parms.ca = ca;
    % lamda
    parms.lamda = lamda;
    
    % constantes de las integrales
    % constante flujo externo
    parms.rkextf = 0;
    % constante single layer
    parms.rksl = 2/(lamda + 1);
    % constante double layer
    parms.rkdl = 2*(lamda - 1)/(lamda + 1);
    
    % constantes de delta de fuerza
    % constante curvatura
    if ka == 1
        parms.rkcurv = 1/ca;
    else
        parms.rkcurv = 0;
    end
    
    % constante gravedad
    if kb == 1    
        parms.rkgrav = 1;
    else
        parms.rkgrav = 0;        
    end
    
    parms.rkelect = 0;
   
    warning('Adimensionalizacion de Bazhlekov gravity termino de gravedad rkgrav = 1 por defecto y  curvatura = 1/bond (definido como capilar)');    
elseif adim == 5 
    % adimensionalizacion de ascoli - gravity
    parms.flow = flow;
    parms.w = 0;
    % parametros
    % capilar
    parms.ca = ca;
    % lamda
    parms.lamda = lamda;
    
    % constantes de las integrales
    % constante flujo externo
    parms.rkextf = 0;
    % constante single layer
    parms.rksl = 2/((lamda + 1)*ca);
    % constante double layer
    parms.rkdl = 2*(lamda - 1)/(lamda + 1);
    
    % constantes de delta de fuerza
    % constante curvatura
    if ka == 1
        parms.rkcurv = 1;
    else
        parms.rkcurv = 0;
    end
    
    % constante gravedad
    if kb == 1    
        parms.rkgrav = ca*(3*(1 + 1.5*lamda)/(1 + lamda));
    else
        parms.rkgrav = 0;        
    end
    
    parms.rkelect = 0;
        
    warning('Adimensionalizacion de Ascoli gravity gravedad rkgrav = ca*(3*(1 + 1.5*lamda)/(1 + lamda))');    
    
end
