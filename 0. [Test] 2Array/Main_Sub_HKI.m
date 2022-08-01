function [InvMat_HKI,VF_For_FFR,p0,Time_Main_Sub_HKI] = Main_Sub_HKI(totfreq,u0,c0,rho0,radius_a,length_l,NN)
%%
tic

%%
VF_nonDim = zeros(length(totfreq),5); VF_Dim = zeros(length(totfreq),5);
InvMat_HKI = cell(length(totfreq),1); VF_For_FFR = cell(length(totfreq),1);

handler = waitbar(0,'Initializing waitbar...'); % waitbar를 띄웁니다.
for NumFreq = 1:length(totfreq)
    freq = totfreq(NumFreq);
    waitbar(NumFreq/length(totfreq),handler,sprintf('Computing... %d Hz', freq));
        
     [GD, GW, p0] = HKI_Sub_SurPres_Mat(freq,u0,c0,rho0,radius_a,length_l,NN);
    
     [InvMat_HKI{NumFreq,1}] = HKI_Sub_CalSurPres(GD, GW);
        
     [VF_For_FFR{NumFreq,1}] = HKI_Sub_CalRadImp(InvMat_HKI{NumFreq,1},radius_a,length_l,NN);

end
close(handler) % 루프가 끝나면 waitbar를 종료한다.

Time_Main_Sub_HKI = toc;
Hour = fix(Time_Main_Sub_HKI/3600); Min = fix(rem(Time_Main_Sub_HKI,3600)/60); Sec = round(rem(rem(Time_Main_Sub_HKI,3600),60));
fprintf('Main_Sub_HKI 계산 소요 시간은 %d시간 %d분 %d초 입니다.\n',Hour,Min,Sec)
end