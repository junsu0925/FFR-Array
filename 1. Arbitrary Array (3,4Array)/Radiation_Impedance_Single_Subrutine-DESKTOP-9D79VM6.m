function [z_r_For_FFR]=Radiation_Impedance_Single_Subrutine(rho0,c0,totfreq,radius_a,length_l,NN)

% Construct Radiation Impednace Matrix
z_r_For_FFR = cell(1,length(totfreq)); % radiation impedance matrix

InvMat_HKI = cell(1,length(totfreq)); 

handler = waitbar(0,'Initializing waitbar...'); % waitbar를 띄웁니다.
for NumFreq = 1:length(totfreq)
    freq = totfreq(NumFreq);
    waitbar(NumFreq/length(freq),handler,sprintf('Computing... %d Hz', totfreq));
        
     [GD, GW] = HKI_Sub_SurPres_Mat(totfreq,c0,rho0,radius_a,length_l,NN);
    
     [InvMat_HKI{1,NumFreq}] = HKI_Sub_CalSurPres(GD, GW);
        
     [z_r_For_FFR{1,NumFreq}] = HKI_Sub_CalRadImp(InvMat_HKI{1,NumFreq},radius_a,length_l,NN);

end
close(handler) % 루프가 끝나면 waitbar를 종료한다.

Time_Main_Sub_HKI = toc;
Hour = fix(Time_Main_Sub_HKI/3600); Min = fix(rem(Time_Main_Sub_HKI,3600)/60); Sec = round(rem(rem(Time_Main_Sub_HKI,3600),60));
fprintf('Main_Sub_HKI 계산 소요 시간은 %d시간 %d분 %d초 입니다.\n',Hour,Min,Sec)

end