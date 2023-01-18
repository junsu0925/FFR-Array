function [z_r_For_FFR,Time_Main_Sub_HKI]=Radiation_Impedance(rho0,c0,totfreq,Geometry,NumR)
tic

Line_Sec = Geometry.Line_Sec;
RealAcuteAngle = Geometry.RealAcuteAngle;
MagDVector = Geometry.MagDVector;

% Construct Radiation Impednace Matrix
z_r_For_FFR = cell(1,length(totfreq)); % radiation impedance matrix

InvMat_HKI = cell(1,length(totfreq)); 
NN = cell(length(totfreq),1);

etaOb = 0.0; % CHIF Point [0.7 0.7 0.7 0.7 0.7 0.7];
ztaOb = 0.0; % CHIF Point [0.0 0.3 -0.3 0.1 -0.1 0.4];

handler = waitbar(0,'Initializing waitbar...'); % waitbar를 띄웁니다.
for NumFreq = 1:length(totfreq)
    freq = totfreq(NumFreq);
    [NN{NumFreq,1}] = HKI_Sub_GenNodeNumber_Ls(MagDVector,freq,c0);
    waitbar(NumFreq/length(totfreq),handler,sprintf('Computing... %d Hz', freq));
        
%      [GD, GW] = HKI_Sub_SurPres_Mat(freq,c0,rho0,radius_a,length_l,NN{NumFreq,1});

%   [GD,GW] = HKI_Sub_SurPres_Mat_Ls(freq,RealAcuteAngle,c0,rho0,Line_Sec,NN{NumFreq,1});

[GD,GW] = HKI_Sub_SurPres_Mat_Ls_CHIEF(freq,RealAcuteAngle,c0,rho0,Line_Sec,NN{NumFreq,1},etaOb,ztaOb);

     [InvMat_HKI{1,NumFreq}] = HKI_Sub_CalSurPres(GD, GW,etaOb,ztaOb);
        
     [z_r_For_FFR{1,NumFreq}] = HKI_Sub_CalRadImp(InvMat_HKI{1,NumFreq},Geometry,NN{NumFreq,1},NumR);

end
close(handler) % 루프가 끝나면 waitbar를 종료한다.

Time_Main_Sub_HKI = toc;
Hour = fix(Time_Main_Sub_HKI/3600); Min = fix(rem(Time_Main_Sub_HKI,3600)/60); Sec = round(rem(rem(Time_Main_Sub_HKI,3600),60));
fprintf('Main_Sub_HKI 계산 소요 시간은 %d시간 %d분 %d초 입니다.\n',Hour,Min,Sec)

end