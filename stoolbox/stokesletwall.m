% calcula la correccion del stokeslet para flujo semiinfinito

function gwc = stokesletwall(nodes,x,w)

numnodes = size(x,1);

% memory assigment
gwc = zeros(3,3,numnodes);

% xpole
nodex0 = nodes;

nodex0_im = [nodex0(1), nodex0(2), 2*w - nodex0(3)];

% green function for wall correction
r = x - repmat(nodex0_im,[numnodes 1]);
rn = normesp(r);
rnq = rn.^3;
rnquin = rn.^5;
r11 = r(:,1).^2;
r12 = r(:,1).*r(:,2);
r13 = r(:,1).*r(:,3);
r22 = r(:,2).^2;
r23 = r(:,2).*r(:,3);
r33 = r(:,3).^2;

gwc(1,1,:) =  (-1./rn - r11./(rnq)) + (2.*nodex0(:,3).*(nodex0(:,3) - r(:,3)).*(1./(rnq) - 3.*r11./rnquin));
gwc(1,2,:) =  (-r12./(rnq)) + (2.*nodex0(:,3).*(nodex0(:,3) - r(:,3)).*(-3.*r12./rnquin));
gwc(1,3,:) =  (-r13./(rnq)) - (2.*nodex0(:,3).*(nodex0(:,3) - r(:,3)).*(-3.*r13./rnquin) - 2.*nodex0(:,3).*(r(:,1)./(rnq)));
gwc(2,1,:) =  (-r12./(rnq)) + (2.*nodex0(:,3).*(nodex0(:,3) - r(:,3)).*(-3.*r12./rnquin));
gwc(2,2,:) =  (-1./rn - r22./(rnq)) + (2.*nodex0(:,3).*(nodex0(:,3) - r(:,3)).*(1./(rnq) -3.*r22./rnquin));
gwc(2,3,:) =  (-r23./(rnq)) - (2.*nodex0(:,3).*(nodex0(:,3) - r(:,3)).*(-3.*r23./rnquin) - 2.*nodex0(:,3).*(r(:,2)./(rnq)));
gwc(3,1,:) =  (-r13./(rnq)) + (2.*nodex0(:,3).*(nodex0(:,3) - r(:,3)).*(-3.*r13./rnquin) - 2.*nodex0(:,3).*(-r(:,1)./(rnq)));
gwc(3,2,:) =  (-r23./(rnq)) + (2.*nodex0(:,3).*(nodex0(:,3) - r(:,3)).*(-3.*r23./rnquin) - 2.*nodex0(:,3).*(-r(:,2)./(rnq)));
gwc(3,3,:) =  (-1./rn - r33./(rnq)) - (2.*nodex0(:,3).*(nodex0(:,3) - r(:,3)).*(1./(rnq) - 3.*r33./rnquin));

% multiplique todo por (-1/8pi)
gwc = (-1/(8*pi)).*gwc;
