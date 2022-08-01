function [ExportObPres,Time_HKI_Sub_CalObPres] = HKI_Sub_CalObPres(totfreq,ro,zo,c0,rho0,radius_a,length_l,NN,p0,InvMat_HKI,UMatrix_for_HKI)
%%
tic

%%
OutP = zeros(length(totfreq),1);
TVR = zeros(length(totfreq),1);

handler = waitbar(0,'Initializing waitbar...'); % waitbar를 띄웁니다.
for NumFreq = 1:length(totfreq)
    freq = totfreq(NumFreq);
    waitbar(NumFreq/length(totfreq),handler,sprintf('Computing... %d Hz', freq));
    
    T_3B_1 = zeros(NN(1),3);
    T_3B_1(:,1) = 1;
    T_3B_2 = zeros(NN(2),3);
    T_3B_2(:,2) = 1;
    T_3B_3 = zeros(NN(3),3);
    T_3B_3(:,3) = 1;
    T_3B = [T_3B_1; T_3B_2; T_3B_3];

    Vu = (T_3B*UMatrix_for_HKI{1,NumFreq});
    [GD_ob, GW_ob] = HKI_Sub_ObPres_Mat(ro,zo,freq,c0,rho0,radius_a,length_l,NN);
    
    OutP(NumFreq,1) = ((GD_ob + GW_ob*InvMat_HKI{NumFreq,1})*Vu);
end
close(handler) % 루프가 끝나면 waitbar를 종료한다.

for NumFreq=1:length(totfreq);
    TVR(NumFreq,1)=20*log10(sqrt(ro^2+zo^2)*sqrt(OutP(NumFreq)*conj(OutP(NumFreq))/2)/1e-6);
end

ExportObPres = [totfreq', OutP, TVR];

Time_HKI_Sub_CalObPres = toc;

Hour = fix(Time_HKI_Sub_CalObPres/3600); Min = fix(rem(Time_HKI_Sub_CalObPres,3600)/60); Sec = round(rem(rem(Time_HKI_Sub_CalObPres,3600),60));
fprintf('HKI_Sub_CalObPres 계산 소요 시간은 %d시간 %d분 %d초 입니다.\n',Hour,Min,Sec)
end