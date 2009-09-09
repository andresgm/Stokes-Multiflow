%% Calcula las derivadas del campo vectorial F(Zita1,Zita2) en la base
%% isoparametrica (Zita1,Zita2), evaluadas en las coordenadas ZitaVect
%% para un elemento triangular plano lineal (3 nodos)
%%      Entradas:
%%          P1,P2,P3 = Componentes del campo vectorial en los nodos
%%                  dim(3,1)
%%              Ex. Pi(3,1) = Valor de F en la comp i en el Nodo 3 LOCAL 
%%          ZitaVect = Arrays de coordenadas isoparametricas de las
%%          posiciones sobre las cuales calcular el campo F.
%%                  dim(2,NumPoints)
%%      Salidas:
%%          dF1Z1,dF1Z2,dF2Z1,dF2Z2,dF3Z1,dF3Z2 = Derivadas dF/dz1, dF/dz2
%%          evaluadas en cada punto de ZitaVect.
%%                 dim(1,NumPoints) 

%% La formulacion parte de F=Phi*Fp donde 
%% Phi = g(Zita1,Zita2) Matriz de funciones de forma y 
%% Fp: Vector de parametros -valores del campo vectorial en los nodos del
%% elemento. Brebbia -Boundary Elements - An introductory course pg 124y175
%% 
%% El valor de df/zi se determina derivando f=phi(zita1,zita2)*fp.

function [df1z1,df1z2,df2z1,df2z2,df3z1,df3z2,df1z3,df2z3,df3z3] = dfzitalin(p1,p2,p3)

    % campo escalar
df1z1 = p1(1,:) - p1(3,:);
df1z2 = p1(2,:) - p1(3,:);
% adicionado derivada respecto zita 3
df1z3 = p1(3,:);
if nargin == 1
df2z1 = [];
df2z2 = [];
df3z1 = [];
df3z2 = [];
% adicionado derivadas respecto zita 3
df2z3 = [];
df3z3 = [];
else
df2z1 = p2(1,:) - p2(3,:);
df2z2 = p2(2,:) - p2(3,:);
df3z1 = p3(1,:) - p3(3,:);
df3z2 = p3(2,:) - p3(3,:);
% adicionado derivadas respecto zita 3
df2z3 = p2(3,:);
df3z3 = p3(3,:);
end
