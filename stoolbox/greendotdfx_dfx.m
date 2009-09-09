%% Calculo del integrando del single layer potential para tension
%% isotropica
% dim(integrandns)= numQudpoints*numElements x 3
% Cada numElements a lo largo de filas corresponden a un k-esimo punto de
% integracion:
% 
% TODO: reorganizar esta rutina

function integrandns = greendotdfx_dfx(greenarrayfcn,deltafatxpole,normalatxpoles,iopt,deltafclose)

numxpoles = size(greenarrayfcn,1);


% [m,n,o] = size(normalatintptVR);


integrandns = cell(size(greenarrayfcn,1),1);



if iopt == 0
    normalatintptl = normalatxpoles;
     % arreglo de valores de normalatintptl y deltafatintpt para ejecutar producto en paralelo
    normalatintptl = repmat(permute(normalatintptl',[1 3 2]),[1 3 1]);
    deltafatxpolepar= repmat(permute(deltafatxpole,[3 2 1]),[3 3 1]);
    % flujo infinito
    % parfor i = 1:size(greenarrayfcn,1)
    parfor i = 1:numxpoles
        % Extraiga las funciones de green del ipole
        array = greenarrayfcn{i};
        % Extraiga el deltaf del ipole
        deltafx0 = deltafatxpole(i);
        % Calcule el deltafx - deltafx0
        deltafx_deltafx0 = deltafatxpolepar - deltafx0;
        % Multiplique deltafx_deltafx0.*normalatintptl
        normdeltf = normalatintptl.*deltafx_deltafx0;
        
        integrandns{i} = sum(array.*normdeltf,1);
    end
elseif iopt == 1
    normalatintptl = normalatxpoles;
     % arreglo de valores de normalatintptl y deltafatintpt para ejecutar producto en paralelo
    normalatintptl = repmat(permute(normalatintptl',[1 3 2]),[1 3 1]);
    % flujo semiinfinito
    parfor i = 1:numxpoles
        % Extraiga las funciones de green del ipole
        array = greenarrayfcn{i};
        % Calcule el deltafx - deltafx0
        deltafx_deltafx0 = deltafatxpole - deltafclose(i);
        % Multiplique deltafx_deltafx0.*normalatintptl
        normdeltf = normalatintptl.*repmat(permute(deltafx_deltafx0,[3 2 1]),[3 3 1]);
        integrandns{i} = sum(array.*normdeltf,1);
    end  
elseif iopt == 2
    normalatintptl = normalatxpoles;
     % arreglo de valores de normalatintptl y deltafatintpt para ejecutar producto en paralelo
    normalatintptl = repmat(permute(normalatintptl',[1 3 2]),[1 3 1]);
    deltafatxpolepar= repmat(permute(deltafatxpole,[3 2 1]),[3 3 1]);
    % integral sin desingularizacion (para marangoni stress)
    % normalatintptl: vector unitario de los esfuerzos de marangoni
    % deltafatxpolepar: magnitud de los esfuerzos de marangoni
    % parfor i = 1:size(greenarrayfcn,1)
    parfor i = 1:numxpoles
        % Extraiga las funciones de green del ipole
        array = greenarrayfcn{i};
        % Multiplique deltafatxpolepar.*normalatintptl
        normdeltf = normalatintptl.*deltafatxpolepar;
        integrandns{i} = sum(array.*normdeltf,1);
    end
elseif iopt == 3
% arreglo de valores de normalatintptl y deltafatintpt para ejecutar producto en paralelo
deltafatxpolel = repmat(permute(deltafatxpole',[1 3 2]),[1 3 1]);
    for i = 1:numxpoles
        % Extraiga las funciones de green del ipole
        array = greenarrayfcn{i};
        % Extraiga el deltafclose
        deltafx0 = deltafclose(i,:);
        % Calcule el deltafx - deltafx0
        deltafx_deltafx0 = deltafatxpolel - repmat(deltafx0',[1 3 numxpoles]);
        integrandns{i} = sum(array.*deltafx_deltafx0,1);
    end
end
