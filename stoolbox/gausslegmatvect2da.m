% Calcula la integral polar 2D de la forma INT(gij*fi)dS usando gauss
% legendre en un elemento triangular plano singular con Origen en el P3 
% del elemento. Realiza una integral doble en dos integrales 1D: 
% respecto a rho(radio) y respecto a phi(Angulo)
%
%   entradas:
%   functionname: @function a evaluar gij
%   vectfld: campo vectorial f definido en cada nodo [F1; F2; F3] 
%   dim(3nodos x 3 comp)
%   triandata.nodes = [P1;P2;P3] dim(3nodos x 3 coord) (polo en P3)
%   triandata.g = metrica del elemento triangular (ver metrictrans)
%   polarparms.xx: parametros de integracion polar de gauss legendre 2D
%   ver gausslegabsweights, gausslegintpt
%
%   Salidas: intval valor de la integral dim(1 x 3)
%
% Tomado de pozrikidis C. pozrikidis, numerical Computation in Science and 
% engineering Chap numerical integration.
% Adaptado de Stokesflow BemliB archivo bump_3d_slp.f rutina
% intgr_trgl_sing y intr_lin_sing

% DESCRIP rutina vectorizada de gauss legendre revizar la vectorizacion

function intval = gausslegmatvect2da(functionname,vectfld,triandata,polarparms)

if nargin < 4
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
  
x1t = triandata.nodo1loc;
x2t = triandata.nodo2loc;
x3t = triandata.nodo3loc;

x1t = permute(x1t,[3 2 1]);
x2t = permute(x2t,[3 2 1]);
x3t = permute(x3t,[3 2 1]);

numptrho=size(wwrho,2);
numptphi=size(ww,2);
numele = size(x1t,3);

piq = 0.25*pi;

% Determine coordenadas punto de integracion en base global x1,x2,x3
tempx1t = repmat(x1t,[numptrho*numptphi 1 1]);
tempx2t = repmat(x2t,[numptrho*numptphi 1 1]);
tempx3t = repmat(x3t,[numptrho*numptphi 1 1]);

x = tempx1t.*repmat(xin,[1 1 numele]) + tempx2t.*repmat(etn,[1 1 numele])...
    + tempx3t.*repmat(ztn,[1 1 numele]);
xlist = reshape(permute(x,[2 1 3]),[3,numptphi*numptrho*numele])';

% interpole los valores del campo vectorial a los puntos de integracion
vectfld1t = repmat(permute(vectfld.fi1loc,[3 2 1]),[numptrho*numptphi 1 1]);
vectfld2t = repmat(permute(vectfld.fi2loc,[3 2 1]),[numptrho*numptphi 1 1]);
vectfld3t = repmat(permute(vectfld.fi3loc,[3 2 1]),[numptrho*numptphi 1 1]);

finterp = vectfld1t.*repmat(xin,[1 1 numele])...
    + vectfld2t.*repmat(etn,[1 1 numele])...
    + vectfld3t.*repmat(ztn,[1 1 numele]);
finterplist = repmat(permute(reshape(permute(finterp,[2 1 3]),[3,numptphi*numptrho*numele]),[1 3 2]),[1 3 1]);

% invoque rutina de funciones de greenwall con polo en x3t

funcionval = feval(functionname,x3t(:,:,1),xlist,0);

funcionval = funcionval{1};

% realice el producto gij*fi por el vector finterp
funcionval = sum(funcionval.*finterplist,1);

% SOlO vAliDO Si numptrho = numptphi
cf = r.*wwrho;                
cftemp1 = permute(reshape(cf,numptrho*numptphi,1),[3 2 1]);
cftemp2 = repmat(cftemp1,[1 3 1]);       
cftemp2 = repmat(cftemp2,[1 1 numele]);

% multiplique las funciones de green por sus pesos de rho
funcionval = funcionval.*cftemp2;

% Sume las las funciones de green para cada phi a lo largo de rho
% multipliquelas por el peso y rmaxh para cada phi y sume
kte_rww = rmaxh.*ww;
kte_rww1 = repmat(repmat(permute(reshape(repmat(kte_rww,[numptphi 1]),[numptphi*numptrho 1]),[3 2 1]),[1 1 numele]),[1 3 1]);
jacarray = repmat(permute(reshape(repmat(piq.*triandata.g',[numptrho*numptphi 1]),[numptrho*numptphi*numele 1]),[3 2 1]),[1 3 1]);
kte_rww1 = kte_rww1.*jacarray;
intval = sum(funcionval.*kte_rww1,3);
            