function [Int_W] = HKI_Sub_IntFn_W(NEq,R,ka,eta_o,eta,theta,zta)
if NEq == 1 % Surface Bottom
    Int_W=eta.*(1+1i*ka.*R(eta,theta))./((R(eta,theta).^3).*exp(1i*ka.*R(eta,theta)));
elseif NEq == 2 % Surface ring
    Int_W=(eta_o.*cos(theta)-1).*(1+1i*ka.*R(zta,theta))./((R(zta,theta).^3).*exp(1i*ka.*R(zta,theta)));
elseif NEq == 3 % Surface Top
    Int_W=eta.*(1+1i*ka.*R(eta,theta))./((R(eta,theta).^3).*exp(1i*ka.*R(eta,theta)));
end
end