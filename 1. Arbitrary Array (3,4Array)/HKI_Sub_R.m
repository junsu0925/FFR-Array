%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate HKI for Axis-Symmetric Cylinder,
% Result : Surface pressure & RadiationImpedance
% : HKI_Sub_R.m
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

function [R]=HKI_Sub_R(eta_o,zta_o,eta,theta,zta)
R = sqrt(eta_o^2 + eta.^2 - 2*eta_o.*eta.*cos(theta) + (zta_o-zta).^2);
end