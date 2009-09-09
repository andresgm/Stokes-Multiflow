% Calcula la integral polar 2D INT(gij)dS usando gauss legendre en un elemento triangular
% plano singular con Origen en el P3 del elemento
% realiza una integral doble en dos integrales 1D: respecto a rho(radio) y
% respecto a phi(Angulo)
%
%   entradas:
%   functionname: @function a evaluar gij
%   triandata.nodes = [P1 P2 P3] dim(3nodos x 3 coord) (polo en P3)
%   triandata.g = metrica del elemento triangular (ver metrictrans)
%   polarparms.xx: parametros de integracion polar de gauss legendre 2D
%   ver gausslegabsweights, gausslegintpt
%
%   Salidas: intval valor de la integral dim(3 x 3)
%
% Tomado de pozrikidis C. pozrikidis, numerical Computation in Science and 
% engineering Chap numerical integration.
% Adaptado de Stokesflow BemliB archivo bump_3d_slp.f rutina
% intgr_trgl_sing y intr_lin_sing


function intval = gausslegmat2d(functionname,triandata,polarparms)

if nargin < 3
    % calcule los parametros de integracion polar 
    % abscisas y pesos de los puntos de integracion gauss-Leg 1D
    % numero de puntos rho - phi 4 (default)
    n = 4;
    [zz,ww] = gausslegabsweights(n);
      % coordenadas y pesos de los puntos de integracion gauss-Leg 2D
    [rmaxh,wwrho,r,xin,etn,ztn] = gausslegintpt(zz,ww);
else
    rmaxh = polarparms.rmaxh;
    wwrho = polarparms.wwrho;
    ww = polarparms.ww;
    xin = polarparms.xin;
    etn = polarparms.etn;
    ztn = polarparms.ztn;
    r = polarparms.r;
end
  
x1t = triandata.nodes(1,:);
x2t = triandata.nodes(2,:);
x3t = triandata.nodes(3,:);
g = triandata.g;

numptrho = size(wwrho,2);
numptphi = size(ww,2);

piq = 0.25*pi;

% Determine coordenadas punto de integracion en base global x1,x2,x3
tempx1t = repmat(x1t,numptrho*numptphi,1);
tempx2t = repmat(x2t,numptrho*numptphi,1);
tempx3t = repmat(x3t,numptrho*numptphi,1);

x = tempx1t.*xin + tempx2t.*etn + tempx3t.*ztn;
% invoque rutina de funciones de greenwall con polo en x3t
funcionval = feval(functionname,x3t,x,0);

funcionval = funcionval{1};

% SOlO vAliDO Si numptrho = numptphi
cf = r.*wwrho;                
cftemp1 = permute(reshape(cf,numptrho*numptphi,1),[3 2 1]);
cftemp2 = repmat(cftemp1,[3 3 1]);                

% multiplique las funciones de green por sus pesos de rho
funcionval = funcionval.*cftemp2;

% Sume las las funciones de green para cada phi a lo largo de rho
% multipliquelas por el peso y rmaxh para cada phi y sume
integralext = zeros(3,3);
minlim = 1;
maxlim = 0;
kte_rww = rmaxh.*ww;
for sd = 1:numptphi
   maxlim = maxlim + numptrho;
   % extraiga funciones de green para un phi dado
   funcionforrho = funcionval(:,:,minlim:1:maxlim);
   % sume las funciones de green para el phi dado 
   integralext = integralext + ...
       sum(funcionforrho.*(kte_rww(sd)),3);
   minlim = maxlim+1;
end

% Completar la integral constante exterior de gauss legendre piq = pi/4 
% y metrica de transformacion g
cf = piq*g;
% integral total  del subelemento singular
intval = integralext.*cf;
            