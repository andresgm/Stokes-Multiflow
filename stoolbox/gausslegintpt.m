% calcula los puntos de integracion en la base xi - eta a partir de los
% puntos en coordenadas polares sobre un triangulo
function [rmaxh,wwrho,r,xin,etn,ztn] = gausslegintpt(zz,ww)
% parametros de gauss legendre
           % Constantes y asignacion de espacio
            piq = 0.25*pi;
    %         pi4 = 4*pi;
%             asm = 0;
        % parametros respecto de phi
        ph = piq.*(1.0 + zz);
        cph = cos(ph);
        sph = sin(ph);
        % r(phi)=1/(cos(phi)+sin(phi))
        rmax = 1./(cph + sph);
        % rmaxh = r(phi)/2
        rmaxh = 0.5.*rmax;
        % parametros respecto de rho
        % Conjunto de posiciones radiales de los puntos de integracion a
        numptphi = size(zz,2);
        numptrho = size(zz,2);
         % Coeficientes peso w.r.t rho
        wwrho = repmat(ww',1,numptphi);   
        temp1 = repmat(rmaxh,numptrho,1);
        temp2 = repmat((1 + zz)',1,numptphi);
        % r=dim(numtptrho,numptphi); radio desde el polo al intpt
        r = temp1.*temp2;
        
        % Determine la posiciones en la base local Zita1 Zita2 de los intpt
        temp1 = repmat(cph,numptrho,1);
		temp2 = repmat(sph,numptrho,1);
        % xi=r*cos(phi)
        xi = r.*temp1;
        % et=r*sin(phi)
        et = r.*temp2;
        % zt=1-Zita1-Zita2
        zt = 1 - xi - et;
        xin = repmat(reshape(xi,numptrho*numptphi,1),1,3);
        etn = repmat(reshape(et,numptrho*numptphi,1),1,3); 
        ztn = repmat(reshape(zt,numptrho*numptphi,1),1,3);