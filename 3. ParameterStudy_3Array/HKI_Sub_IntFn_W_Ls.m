%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate HKI for Axis-Symmetric Cylinder,
% Result : Surface pressure & RadiationImpedance
% : HKI_Sub_IntFn_W.m
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

function [Int_W] = HKI_Sub_IntFn_W_Ls(NEq,R,kLmin,psi,Ref_etaztan,eta_o,zta_o,IntRange,NorDVector,Ls,theta)

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
    Int_W=eta.*(((1+1i*kLmin.*R(eta,theta,zta))./...
        ((R(eta,theta,zta).^3).*exp(1i*kLmin.*R(eta,theta,zta)))).*...
        (zta_o-zta).*NorDVector(2));
elseif NEq == 2 % Case 2
    Int_W=eta.*(((1+1i*kLmin.*R(eta,theta,zta))./...
        ((R(eta,theta,zta).^3).*exp(1i*kLmin.*R(eta,theta,zta)))).*...
        (((eta_o-eta.*cos(theta)).*(NorDVector(1).*cos(theta)))+...
        ((-eta.*sin(theta)).*(NorDVector(1).*sin(theta)))));
elseif NEq == 3 % Case 3    
    Int_W=((eta_1+eta_2)./2).*...
        (((1+1i*kLmin.*R(eta,theta,zta))./...
        ((R(eta,theta,zta).^3).*exp(1i*kLmin.*R(eta,theta,zta)))).*...
        (((eta_o-eta.*cos(theta)).*(NorDVector(1).*cos(theta)))+...
        ((-eta.*sin(theta)).*(NorDVector(1).*sin(theta)))+...
        ((zta_o-zta).*NorDVector(2))));
end
end