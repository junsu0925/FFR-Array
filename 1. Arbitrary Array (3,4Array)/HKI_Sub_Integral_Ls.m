%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate HKI for Axis-Symmetric Cylinder,
% Result : Surface pressure & RadiationImpedance
% : HKI_Sub_Integral.m
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

function IntResult = HKI_Sub_Integral_Ls(Int_Fn,IntRange,StartThata,EndTheta,IntError)
Method = 'auto';
% IntResult=integral2(Int_SFn,IntRange(1)+IntError,IntRange(2)-IntError,0,2*pi,'Method',Method,'AbsTol',1e-12,'RelTol',1e-8);
% IntResult=integral2(Int_Fn,IntRange(1)+IntError,IntRange(2)-IntError,0,2*pi,'Method',Method);
IntResult=integral2(Int_Fn,IntRange(1)+IntError,IntRange(2)-IntError,StartThata,EndTheta,'Method',Method);
end