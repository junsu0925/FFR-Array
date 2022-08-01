function [Int_D] = HKI_Sub_IntFn_D(NEq,R,ka,eta,theta,zta)
if NEq == 1 % Surface Bottom
    Int_D=eta.*exp(-1i.*ka.*R(eta,theta))./R(eta,theta);
elseif NEq == 2 % Surface ring
    Int_D=exp(-1i.*ka.*R(zta,theta))./R(zta,theta);
elseif NEq == 3 % Surface Top
    Int_D=eta.*exp(-1i.*ka.*R(eta,theta))./R(eta,theta);
end
end