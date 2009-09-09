%% green's functions and stress Tensors for semi-infinite flow bounded by a
%% plane wall located at x1 = w
%%          Entradas: 
%%              nodes: Collocation Points, Polos o nodos del problema
%%              x: Puntos x de la funci?n de green (normalmente 
%%                  Puntos de integraci?n)         
%%              iopt: 1 calcule stress tensor; n.e 1 calcule solo
%%                  stockeslet
%%              w: distancia de la pared en direcci?n x1 desde la base
%%                  global
%%          salidas:
%%              greenWallfcn: Cell array de funciones de green
%%                  dimCell(numnodes,1)
%%                  Cada k Cell Array contiene numxpoints matrices de 3 x 3
%%              stressTensorWallfcn: Cell array de Tensor de esfuerzos

%% gt(x,x0) = ginf + gwall

%% Tomado de Pozrikidis Boundary integral and singularity methods for
%% Linearized Viscous flow Chap 3 pg 84

% TODO: enlazar stress tensor. Potencialment cambiar el nombre
function [greeninffcn,closenode] = greeninfp(nodes,x,iopt)

numnodes = size(nodes,1);

% memory assigment
greeninffcn = cell(numnodes,1);
closenode = zeros(numnodes,1);

%     parfor i = 1:numnodes
    for i = 1:numnodes
        % xpole
        nodex0 = nodes(i,:);
        % green function for free space part
        [greeninffcn{i}, closenode(i)] = stokesletp(nodex0,x);

        if iopt == 1 % stress Tensor Calculation

        end
    end
