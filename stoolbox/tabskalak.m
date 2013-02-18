function tab = ...
    tabskalak(Aabctr,JAabctr,aab,aabctr,jAab,gssk,csk)
% Returns the contravariant representation of the Cauchy tension tensor
    
    I1 = ...
        (Aabctr(1,:).*aab(1,:) + Aabctr(2,:).*aab(2,:) + ...
        2*Aabctr(3,:).*aab(3,:)-2)';
    I2 = JAabctr.*jAab-1;
    
    % Skalak constitutive equation
    dwdI1 = gssk*(I1+1)/4;
    dwdI2 = (csk*I2-gssk)/4;
    
    Js = sqrt(JAabctr.*jAab);
    
    tab = ...
        2*((repmat(dwdI1,1,3).*Aabctr')./repmat(Js,1,3) + ...
        repmat(Js,1,3).*repmat(dwdI2,1,3).*aabctr');
end