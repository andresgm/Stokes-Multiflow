% calcula los terminos de la evolucion de la concentracion de sulfactantes
% en una superficie.
% Basado en Bazhlekov I.B "Boundary Integral Method for Deformable
% Interfaces in the Presence of Insoluble Surfactants" I. Lirkov et al.
% (Eds.): LSSC 2003, LNCS 2907, pp. 355-362, 2004.

% struct: estructura de la geometria de la malla
% u: Velocidad Hidrodinamica en la interfase
% w: Velocidad tangencial arbitraria (Vel Adapt Mesh)
% gamma: campo escalar de concentracion
% pe: numero de Peclet de la concentracion
% surfopt: opciones de surfactants
function [dgammadt,terms] = surfactexp(struct,gamma,u,w,pe,surfopt)


% preproceso de las entradas
numnodes = size(struct.nodes,1);
numelements = size(struct.elements,1);

 
elefield = isfield(struct,'element2node');
if elefield == 0
    % calcule los elementos vecinos al cada nodo
    struct.element2node = element2node(struct.elements);
end

%% Generales
% Calcule el tensor de proyeccion tangencial en cada nodo
pjnode = projtensor(struct.normal);
% calcule la proyeccion normal del campo de velocidades u
u_nmag = sum(u.*struct.normal,2);
u_n = repmat(u_nmag,[1 3]).*struct.normal;
% calcule el campo de velocidades hidrodinamico tangencial
u_s = u - u_n;

u_s1 = zeros(numnodes,3);
for k=1:numnodes
u_s1(k,:) = (pjnode(:,:,k)*u(k,:)')';
end

    %% Calculo del primer termino (w + u_s).grad_s(gamma)
% invoque rutina de transporte de un campo escalar (surface)

if surfopt.opt == 1 
% los nodos solo se mueven con la velocidad normal
term1 = zeros(numnodes,1);
elseif surfopt.opt == 2 
% los nodos se mueven con la velocidad normal y adaptacion
term1 = transscfld(struct,gamma,w);
elseif surfopt.opt == 3 
% los nodos se mueven con la velocidad hidrodinamica
term1 = transscfld(struct,gamma,u_s);
elseif surfopt.opt == 4 
% los nodos se mueven con la velocidad hidrodinamica y adaptacion
term1 = transscfld(struct,gamma,(w + u_s));
end

%% Calculo del tercer termino gamma.curv<u,n>
term3 = gamma.*struct.curv.*u_nmag;

%% Calculo del cuarto termino lap_s(gamma)
% invoque laplace beltrami de un campo escalar
term4 = lapbelscfld(struct,gamma);

%% Calculo del segundo termino div_s(gamma.u_s)

% calcule gamma.u_s en los nodos
gammau_s = u_s.*repmat(gamma,[1 3]);

% Coordenadas nodo 1, 2 y 3 vertices
n1 = struct.nodes(struct.elements(:,1),:);
n2 = struct.nodes(struct.elements(:,2),:);
n3 = struct.nodes(struct.elements(:,3),:);
% Coordenadas nodos intermedios 4,5,6
n4 = (struct.nodes(struct.elements(:,1),:) + struct.nodes(struct.elements(:,2),:)).*0.5;
n5 = (struct.nodes(struct.elements(:,3),:) + struct.nodes(struct.elements(:,2),:)).*0.5;
n6 = (struct.nodes(struct.elements(:,3),:) + struct.nodes(struct.elements(:,1),:)).*0.5;
% Coordenadas nodo centroidad 7.
n7 = struct.nodes(struct.elements(:,1),:).*(1/3) + ...
    struct.nodes(struct.elements(:,2),:).*(1/3) + ...
    struct.nodes(struct.elements(:,3),:).*(1/3);

% Diferencia entre nodo 7 e intermedios (genericamente n7_ni)
n7_n4 = n7 - n4;
n7_n5 = n7 - n5;
n7_n6 = n7 - n6;

n7_n1 = n7 - n1;
n7_n2 = n7 - n2;
n7_n3 = n7 - n3;

% Normas o longitudes de los segmentos en cada elemento
Sa = normesp(n7_n4);
Sb = normesp(n7_n5);
Sc = normesp(n7_n6);

% gammaus en nodo 4,5 y 6 intermedios
gammau_s4 = (gammau_s(struct.elements(:,1),:) + gammau_s(struct.elements(:,2),:)).*0.5;
gammau_s5 = (gammau_s(struct.elements(:,3),:) + gammau_s(struct.elements(:,2),:)).*0.5;
gammau_s6 = (gammau_s(struct.elements(:,3),:) + gammau_s(struct.elements(:,1),:)).*0.5;
% calcule gamma.u_s en el centroide se los elementos
gammau_s7 = gammau_s(struct.elements(:,1),:).*(1/3) + ...
    gammau_s(struct.elements(:,2),:).*(1/3) + ...
    gammau_s(struct.elements(:,3),:).*(1/3);

% % asigne espacio al termino 2
% term2 = zeros(struct.numnodes,1);
% 
% for jnode = 1:numnodes
%     % extraiga los elementos del jnode
%     eleind = struct.element2node{jnode};
%     % arreglo termino2
%     intl = zeros(size(eleind,2),1);
%     
%     % extraiga los nodos de los elementos vecinos al jnodes
%     indnod = struct.elements(eleind,:);
%     % extraiga las normales a los elementos
%     normalele = struct.normalele(eleind,:);
%         
%     % caso 1 cuando jnode es el nodo 1 local
%     tempind = indnod(:,1) == jnode;
%     if sum(tempind,1) > 0
%         % 3 nodo4
%         s1 = Sa(eleind(tempind));
% 
%         t1 = n7_n4(eleind(tempind),:)./repmat(s1,[1 3]);
% 
%         b1 = cross(permute(t1,[3 2 1]),permute(normalele(tempind,:),[3 2 1]));
%         b1 = permute(b1,[3 2 1]);
%         
%         % comprueba que es vector exterior
%         comprobar = sum(b1.*n7_n1(eleind(tempind),:),2);
%         
%         sss = find(comprobar < 0);
%         if size(sss,1) ~= 0
%         comprobar
%         end    
%         % 1 nodo 6
%         s2 = Sc(eleind(tempind));
% 
%         t2 = -n7_n6(eleind(tempind),:)./repmat(s2,[1 3]);
% 
%         b2 = cross(permute(t2,[3 2 1]),permute(normalele(tempind,:),[3 2 1]));
%         b2 = permute(b2,[3 2 1]);
%         
%         % comprueba que es vector exterior
%         comprobar = sum(b2.*n7_n1(eleind(tempind),:),2);
% 
%         sss = find(comprobar < 0);
%         if  size(sss,1) ~= 0
%         comprobar
%         end    
%         
%         % integre los segmentos del termino 2
%         intl(tempind,:) = sum(b1.*(gammau_s4(tempind,:) + gammau_s7(tempind,:)).*repmat(s1,[1 3]).*0.5,2)...
%             + sum(b2.*(gammau_s7(tempind,:) + gammau_s6(tempind,:)).*repmat(s2,[1 3]).*0.5,2);
%         
%     end
% 
%     % caso 2 cuando en 2 esta el jnode
%     tempind = indnod(:,2) == jnode;
%     if sum(tempind,1) > 0
%         % 1 nodo 4
%         s1 = Sa(eleind(tempind));
% 
%         t1 = -n7_n4(eleind(tempind),:)./repmat(s1,[1 3]);
% 
%         b1 = cross(permute(t1,[3 2 1]),permute(normalele(tempind,:),[3 2 1]));
%         b1 = permute(b1,[3 2 1]);
% 
%         % comprueba que es vector exterior
%         comprobar = sum(b1.*n7_n2(eleind(tempind),:),2);
%         
%         sss = find(comprobar < 0);
%         if size(sss,1) ~= 0
%         comprobar
%         end  
%         
%         % 2 nodo5
%         s2 = Sb(eleind(tempind));
% 
%         t2 = n7_n5(eleind(tempind),:)./repmat(s2,[1 3]);
% 
%         b2 = cross(permute(t2,[3 2 1]),permute(normalele(tempind,:),[3 2 1]));
%         b2 = permute(b2,[3 2 1]);
%         
%         % comprueba que es vector exterior
%         comprobar = sum(b2.*n7_n2(eleind(tempind),:),2);
% 
%         sss = find(comprobar < 0);
%         if size(sss,1) ~= 0
%         comprobar
%         end  
% 
%         % integre los segmentos del termino 2
%         intl(tempind,:) = sum(b1.*(gammau_s4(tempind,:) + gammau_s7(tempind,:)).*repmat(s1,[1 3]).*0.5,2)...
%             + sum(b2.*(gammau_s5(tempind,:) + gammau_s7(tempind,:)).*repmat(s2,[1 3]).*0.5,2);
%     end    
% 
% 
%     % caso 3 cuando en 3 esta el jnode
%     tempind = indnod(:,3) == jnode;
%     if sum(tempind,1) > 0 
%         % 3 nodo6
%         s1 = Sc(eleind(tempind));
% 
%         t1 = n7_n6(eleind(tempind),:)./repmat(s1,[1 3]);
% 
%         b1 = cross(permute(t1,[3 2 1]),permute(normalele(tempind,:),[3 2 1]));
%         b1 = permute(b1,[3 2 1]);
% 
%         % comprueba que es vector exterior
%         comprobar = sum(b1.*n7_n3(eleind(tempind),:),2);
% 
%         sss = find(comprobar < 0);
%         if size(sss,1) ~= 0
%         comprobar
%         end  
%         
%         % 2 nodo5
%         s2 = Sb(eleind(tempind));
% 
%         t2 = -n7_n5(eleind(tempind),:)./repmat(s2,[1 3]);
% 
%         b2 = cross(permute(t2,[3 2 1]),permute(normalele(tempind,:),[3 2 1]));
%         b2 = permute(b2,[3 2 1]);
% 
%         % comprueba que es vector exterior
%         comprobar = sum(b2.*n7_n3(eleind(tempind),:),2);
% 
%         sss = find(comprobar < 0);
%         if size(sss,1) ~= 0
%         comprobar
%         end  
%         
%         % integre los segmentos del termino 2
%         intl(tempind,:) = sum(b1.*(gammau_s6(tempind,:) + gammau_s7(tempind,:)).*repmat(s1,[1 3]).*0.5,2)...
%             + sum(b2.*(gammau_s5(tempind,:) + gammau_s7(tempind,:)).*repmat(s2,[1 3]).*0.5,2);
%     end
% 
%     % calcule el termino 2 resultante
%     term2(jnode) = sum(intl,1)./struct.dsi(jnode);
% end
% 
% %% dgammadt y terminos individuales
% terms = [term1 term2 term3 term4];
% dgammadt = term1 - term2 - 2.*term3 + (1/pe).*term4;


%% Calculo simple del segundo termino

gus = repmat(gamma,[1 3]).*u_s;
int2t = zeros(numnodes,1);
tic
for k=1:numnodes
    % extraiga los elementos del nodo k
    elek = struct.element2node{k};
    
    
    for z =1:size(elek,2)
       % extraiga los nodos del elemento z
       nodz = struct.elements(elek(z),:);
       
       % busque la posicion del nodo k en el elemento
       posz = find(nodz == k);
       
       if posz == 1
           % el nodo esta en el nodo 1 local
           
           nodo4 = (struct.nodes(nodz(1),:) + struct.nodes(nodz(2),:)).*0.5;
           nodo7 = struct.nodes(nodz(1),:).*(1/3) + struct.nodes(nodz(2),:).*(1/3) + ...
                   struct.nodes(nodz(3),:).*(1/3);
           nodo6 = (struct.nodes(nodz(3),:) + struct.nodes(nodz(1),:)).*0.5;
           
           gus4 = (gus(nodz(1),:) + gus(nodz(2),:)).*0.5;
           gus7 = gus(nodz(1),:).*(1/3) + gus(nodz(2),:).*(1/3) + ...
                   gus(nodz(3),:).*(1/3);
           gus6 = (gus(nodz(3),:) + gus(nodz(1),:)).*0.5;
           
           % vector tangencial b1
           n7_n4 = nodo7 - nodo4;
           l1 = NormEsp(n7_n4);
           % vector normal al segmento b1
           b1 = cross(n7_n4,struct.normalele(elek(z),:));
           % normalice el vector b1 (ojo se quita si acaso
           b1 = b1./repmat(NormEsp(b1),[1 3]);
           
           % porcion 1 de la integral del elemento
           int1 = sum(l1.*0.5.*(gus4 + gus7) .* b1,2);
           
                                 
           % vector tangencial b2
           n6_n7 = nodo6 - nodo7;
           l2 = NormEsp(n6_n7);
           % vector normal al segmento b1
           b2 = cross(n6_n7,struct.normalele(elek(z),:));
           % normalice el vector b1 (ojo se quita si acaso
           b2 = b2./repmat(NormEsp(b2),[1 3]);
           
           % porcion 1 de la integral del elemento
           int2 = sum(l2.*0.5.*(gus6 + gus7) .* b2,2);
            
           % integral del elemento
           intele = int1 + int2;
           
       elseif posz == 2
           % el nodo esta en el nodo 2 local
           
           nodo4 = (struct.nodes(nodz(1),:) + struct.nodes(nodz(2),:)).*0.5;
           nodo7 = struct.nodes(nodz(1),:).*(1/3) + struct.nodes(nodz(2),:).*(1/3) + ...
                   struct.nodes(nodz(3),:).*(1/3);
           nodo5 = (struct.nodes(nodz(3),:) + struct.nodes(nodz(2),:)).*0.5;
           
           gus4 = (gus(nodz(1),:) + gus(nodz(2),:)).*0.5;
           gus7 = gus(nodz(1),:).*(1/3) + gus(nodz(2),:).*(1/3) + ...
                   gus(nodz(3),:).*(1/3);
           gus5 = (gus(nodz(3),:) + gus(nodz(2),:)).*0.5;
           
           % vector tangencial b1
           n4_n7 = nodo4 - nodo7;
           l1 = NormEsp(n4_n7);
           % vector normal al segmento b1
           b1 = cross(n4_n7,struct.normalele(elek(z),:));
           % normalice el vector b1 (ojo se quita si acaso
           b1 = b1./repmat(NormEsp(b1),[1 3]);
           
           % porcion 1 de la integral del elemento
           int1 = sum(l1.*0.5.*(gus4 + gus7) .* b1,2);
           
                                 
           % vector tangencial b2
           n7_n5 = nodo7 - nodo5;
           l2 = NormEsp(n7_n5);
           % vector normal al segmento b1
           b2 = cross(n7_n5,struct.normalele(elek(z),:));
           % normalice el vector b1 (ojo se quita si acaso
           b2 = b2./repmat(NormEsp(b2),[1 3]);
           
           % porcion 1 de la integral del elemento
           int2 = sum(l2.*0.5.*(gus5 + gus7) .* b2,2);
            
           % integral del elemento
           intele = int1 + int2;
           
       elseif posz == 3
           % el nodo esta en el nodo 3 local
           
           nodo6 = (struct.nodes(nodz(1),:) + struct.nodes(nodz(3),:)).*0.5;
           nodo7 = struct.nodes(nodz(1),:).*(1/3) + struct.nodes(nodz(2),:).*(1/3) + ...
                   struct.nodes(nodz(3),:).*(1/3);
           nodo5 = (struct.nodes(nodz(3),:) + struct.nodes(nodz(2),:)).*0.5;
           
           gus6 = (gus(nodz(1),:) + gus(nodz(3),:)).*0.5;
           gus7 = gus(nodz(1),:).*(1/3) + gus(nodz(2),:).*(1/3) + ...
                   gus(nodz(3),:).*(1/3);
           gus5 = (gus(nodz(3),:) + gus(nodz(2),:)).*0.5;
           
           % vector tangencial b1
           n7_n6 = nodo7 - nodo6;
           l1 = NormEsp(n7_n6);
           % vector normal al segmento b1
           b1 = cross(n7_n6,struct.normalele(elek(z),:));
           % normalice el vector b1 (ojo se quita si acaso
           b1 = b1./repmat(NormEsp(b1),[1 3]);
           
           % porcion 1 de la integral del elemento
           int1 = sum(l1.*0.5.*(gus6 + gus7) .* b1,2);
           
                                 
           % vector tangencial b2
           n5_n7 = nodo5 - nodo7;
           l2 = NormEsp(n5_n7);
           % vector normal al segmento b1
           b2 = cross(n5_n7,struct.normalele(elek(z),:));
           % normalice el vector b1 (ojo se quita si acaso
           b2 = b2./repmat(NormEsp(b2),[1 3]);
           
           % porcion 1 de la integral del elemento
           int2 = sum(l2.*0.5.*(gus5 + gus7) .* b2,2);
            
           % integral del elemento
           intele = int1 + int2;
           
       end    
       
       int2t(k) = int2t(k) + intele;
    end
    
end
term2 = int2t./struct.dsi;
terms = [term1 term2 term3 term4];
dgammadt = term1 - term2 - 2.*term3 + (1/pe).*term4;

