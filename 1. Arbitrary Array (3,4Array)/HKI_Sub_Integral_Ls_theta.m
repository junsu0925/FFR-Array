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

function IntResult = HKI_Sub_Integral_Ls_theta(Int_Fn,IntRange,IntTheta,IntError)
% Method = 'auto';
% IntResult=integral2(Int_Fn,IntRange(1)+IntError,IntRange(2)-IntError,IntTheta(1),IntTheta(2),'Method',Method);
% IntResult=IntResult+integral2(Int_Fn,IntRange(1)+IntError,IntRange(2)-IntError,IntTheta(2),IntTheta(3),'Method',Method);
% IntResult=IntResult+integral2(Int_Fn,IntRange(1)+IntError,IntRange(2)-IntError,IntTheta(3),IntTheta(4),'Method',Method);
% IntResult=IntResult+integral2(Int_Fn,IntRange(1)+IntError,IntRange(2)-IntError,IntTheta(4),IntTheta(5),'Method',Method);
% IntResult=IntResult+integral2(Int_Fn,IntRange(1)+IntError,IntRange(2)-IntError,IntTheta(5),IntTheta(6),'Method',Method);
% IntResult=IntResult+integral2(Int_Fn,IntRange(1)+IntError,IntRange(2)-IntError,IntTheta(6),IntTheta(7),'Method',Method);
% IntResult=IntResult+integral2(Int_Fn,IntRange(1)+IntError,IntRange(2)-IntError,IntTheta(7),IntTheta(8),'Method',Method);
% IntResult=IntResult+integral2(Int_Fn,IntRange(1)+IntError,IntRange(2)-IntError,IntTheta(8),IntTheta(9),'Method',Method);

% Method = 'auto';
% IntResult=integral2(Int_Fn,IntRange(1)+IntError,IntRange(2)-IntError,IntTheta(1),IntTheta(2),'Method',Method);
% IntResult=IntResult+integral2(Int_Fn,IntRange(1),IntRange(2),IntTheta(2),IntTheta(3),'Method',Method);
% IntResult=IntResult+integral2(Int_Fn,IntRange(1),IntRange(2),IntTheta(3),IntTheta(4),'Method',Method);
% IntResult=IntResult+integral2(Int_Fn,IntRange(1),IntRange(2),IntTheta(4),IntTheta(5),'Method',Method);
% IntResult=IntResult+integral2(Int_Fn,IntRange(1),IntRange(2),IntTheta(5),IntTheta(6),'Method',Method);
% IntResult=IntResult+integral2(Int_Fn,IntRange(1),IntRange(2),IntTheta(6),IntTheta(7),'Method',Method);
% IntResult=IntResult+integral2(Int_Fn,IntRange(1),IntRange(2),IntTheta(7),IntTheta(8),'Method',Method);
% IntResult=IntResult+integral2(Int_Fn,IntRange(1)+IntError,IntRange(2)-IntError,IntTheta(8),IntTheta(9),'Method',Method);

%% Ver2
numDivLs = 2;

CentAngle = 40;
numCenTheta = 20;
numSideTheta = 64;

divCenTheta = CentAngle/numCenTheta;
divSideTheta = (360-CentAngle)/numSideTheta;
Int_Firtheta = deg2rad((0:divCenTheta:CentAngle/2));
Int_Midtheta = deg2rad((CentAngle/2+divSideTheta:divSideTheta:360-CentAngle/2));
Int_Fintheta = deg2rad((360-CentAngle/2+divCenTheta:divCenTheta:360));
Int_theta = [Int_Firtheta Int_Midtheta Int_Fintheta];

divLs = (IntRange(2)- IntRange(1))/numDivLs;
Int_Ls = (IntRange(1):divLs:IntRange(2));
IntResult = 0;

for numj = 1:length(Int_theta)-1
    for numi = 1:length(Int_Ls)-1
        Fn(numi,numj) = Int_Fn((Int_Ls(numi)+Int_Ls(numi+1))/2,(Int_theta(numj)+Int_theta(numj+1))/2);
        DeltaArea = divLs*(Int_theta(numj+1) - Int_theta(numj));
        IntResult = IntResult + Fn(numi,numj)*DeltaArea;
    end
end
%% Ver1
% numDivLs = 7;
% % numTheta = 180;
% 
% % divTheta = 360/numTheta;
% % Int_theta = deg2rad((1:divTheta:360));
% 
% CentAngle = 40;
% numCenTheta = 20;
% numSideTheta = 20;
% 
% divCenTheta = CentAngle/numCenTheta;
% divSideTheta = (360-CentAngle)/numSideTheta;
% Int_Firtheta = deg2rad((0:divCenTheta:CentAngle/2));
% Int_Midtheta = deg2rad((CentAngle/2+divSideTheta:divSideTheta:360-CentAngle/2));
% Int_Fintheta = deg2rad((360-CentAngle/2+divCenTheta:divCenTheta:360));
% Int_theta = [Int_Firtheta Int_Midtheta Int_Fintheta];
% 
% divLs = ((IntRange(2)-IntError) - (IntRange(1)+IntError))/numDivLs;
% Int_Ls = ((IntRange(1)+IntError):divLs:(IntRange(2)-IntError));
% IntResult = 0;
% 
% for numj = 1:length(Int_theta)-1
%     for numi = 1:length(Int_Ls)-1
%         Fn(numi,numj) = Int_Fn((Int_Ls(numi)+Int_Ls(numi+1))/2,(Int_theta(numj)+Int_theta(numj+1))/2);
%         DeltaArea = divLs*(Int_theta(numj+1) - Int_theta(numj));
%         IntResult = IntResult + Fn(numi,numj)*DeltaArea;
%         
% %         Fn1(numi,numj) = Int_Fn(Int_Ls(numi),Int_theta(numj));
% %         Fn2(numi,numj) = Int_Fn(Int_Ls(numi+1),Int_theta(numj));
% %         Fn3(numi,numj) = Int_Fn(Int_Ls(numi),Int_theta(numj+1));
% %         Fn4(numi,numj) = Int_Fn(Int_Ls(numi+1),Int_theta(numj+1));
% %         
% %         IntResult = IntResult + mean([Fn1(numi,numj),Fn2(numi,numj),Fn3(numi,numj),Fn4(numi,numj)])*DeltaArea;
%     end
% end

% IntResult
% IntResult_2

% IntResult=integral2(Int_SFn,IntRange(1)+IntError,IntRange(2)-IntError,0,2*pi,'Method',Method,'AbsTol',1e-12,'RelTol',1e-8);
% IntResult=integral2(Int_Fn,IntRange(1)+IntError,IntRange(2)-IntError,0,2*pi,'Method',Method);
end