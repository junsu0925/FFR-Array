%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate HKI for Axis-Symmetric Cylinder,
% Result : Surface pressure & RadiationImpedance
% : HKI_Sub_IntFn_D.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wonkyu Moon
% Department of Mechanical Engineering
% Pohang University of Science and Technology
% Pohang, Korea
% wkmoon@postech.ac.kr
% and
% Kyounghun Been
% Department of Mechanical Engineering
% Pohang University of Science and Technology
% Pohang, Korea
% khbeen@postech.ac.kr

function [Int_D] = HKI_Sub_IntFn_D_Ls(NEq,R,kLmin,psi,Ref_etaztan,IntRange,NorDVector,Ls,theta)

eta = Ref_etaztan(1,1) + Ls*cos(psi);
zta = Ref_etaztan(1,2) + Ls*sin(psi);

% Ver 4
eta_1 = (Ref_etaztan(1,1) + IntRange(1)*cos(psi))*2;
eta_2 = 2*(Ls-IntRange(1))*cos(psi);

% Ver 3
% eta_1 = (Ref_etaztan(1,1) + IntRange(1)*cos(psi))*2;
% eta_2 = 2*Ls*cos(psi);

% Ver 2
% eta_1 = Ref_etaztan(1,1) + IntRange(1)*cos(psi);
% eta_2 = Ref_etaztan(1,1) + Ls*cos(psi);

% Ver 1
% eta_1 = Ref_etaztan(1,1) + IntRange(1)*cos(psi);
% eta_2 = Ref_etaztan(1,1) + IntRange(2)*cos(psi);

if NEq == 1 % Case 1
    Int_D=eta.*(exp(-1i.*kLmin.*R(eta,theta,zta))./R(eta,theta,zta)).*(NorDVector(2))^2;
%     Int_D=eta.*(exp(-1i.*kLmin.*R(eta,theta,zta))./R(eta,theta,zta)).*NorDVector(2);
elseif NEq == 2 % Case 2
    Int_D=eta.*(exp(-1i.*kLmin.*R(eta,theta,zta))./R(eta,theta,zta)).*(NorDVector(1))^2;
%     Int_D=eta.*(exp(-1i.*kLmin.*R(eta,theta,zta))./R(eta,theta,zta)).*NorDVector(1);
elseif NEq == 3 % Case 3    
    Int_D=((eta_1+eta_2)./2).*...
        (exp(-1i.*kLmin.*R(eta,theta,zta))./R(eta,theta,zta)).*...
        ((NorDVector(1))*NorDVector(1)+(NorDVector(2))*NorDVector(2));
end
end

% function [Int_D] = HKI_Sub_IntFn_D(NEq,R,ka,eta,theta,zta)
% if NEq == 1 % Surface Bottom
%     Int_D=eta.*exp(-1i.*ka.*R(eta,theta))./R(eta,theta);
% elseif NEq == 2 % Surface ring
%     Int_D=exp(-1i.*ka.*R(zta,theta))./R(zta,theta);
% elseif NEq == 3 % Surface Top
%     Int_D=eta.*exp(-1i.*ka.*R(eta,theta))./R(eta,theta);
% end
% end