%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate HKI for Axis-Symmetric Cylinder,
% Result : Surface pressure & RadiationImpedance
% : HKI_Sub_SFn.m
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

function [SFn]=HKI_Sub_SFn(SFnNum,NodeRange,eta,IntRange)
if IntRange(1) < min(NodeRange) || IntRange(2) > max(NodeRange)
    SFn = 0;
elseif SFnNum == 1
    SFn = (eta-IntRange(1))/(IntRange(2)-IntRange(1)); % »ó½Â, x(num-1) = a < x < x(num) = b
elseif SFnNum == 2
    SFn = (IntRange(2)-eta)/(IntRange(2)-IntRange(1)); % ÇÏ°­, x(num) = b < x < x(num+1) = c
end
end