%%
clc
clear
close all

%% MultiThread 설정
% NumOfCore = 6;
% parpool(NumOfCore);

%%
Nb = 21; Nr = 35; Ng = 21; Nt = 21;
NN = [Nb, Nr, Ng, Nt];

% freq = 1000:100:4000; % Calculation Frequency
freq = 1000; % Calculation Frequency
% freq = 1000:100:4000; % Calculation Frequency
u0 = 1; % Define initial Velocity

c_water = 1500;
rho_water = 999;

radius_a = 0.175;% Ring Average Radius
length_l = 0.1925; % Ring Height


%%
ro = 1; % the r-distance to calculate TVR
zo = 0; % the z-distance to calculate TVR

%%
[InvMat_HKI,VF_For_FFR,p0,Time_Main_Sub_HKI] = ...
    Main_Sub_HKI(freq,u0,c_water,rho_water,radius_a,length_l,NN);

[ExportCircuit,UMatrix_for_HKI,Time_Main_Sub_Single_FFR] =...
    Main_Sub_Single_FFR(freq,c_water,rho_water,radius_a,length_l,VF_For_FFR);

% [ExportObPres,Time_HKI_Sub_CalObPres] =...
%     HKI_Sub_CalObPres(freq,ro,zo,c_water,rho_water,radius_a,length_l,NN,p0,InvMat_HKI,UMatrix_for_HKI);

%%
Total_Time = Time_Main_Sub_HKI + Time_Main_Sub_Single_FFR + Time_HKI_Sub_CalObPres;
Hour = fix(Total_Time/3600); Min = fix(rem(Total_Time,3600)/60); Sec = round(rem(rem(Total_Time,3600),60));
fprintf('전체 계산 소요 시간은 %d시간 %d분 %d초 입니다.\n',Hour,Min,Sec)

%%
figure(1)
plot(ExportCircuit(:,2),real(ExportCircuit(:,3)),'LineWidth',2)
grid on
xlabel('ka/ka_0','fontsize',20, 'fontangle','italic');
ylabel('Conductance','fontsize',20, 'fontangle','italic');
set(gca, 'fontsize',16)
set(gcf, 'color', 'w')

figure(2)
plot(ExportCircuit(:,2),ExportObPres(:,3),'LineWidth',2)
grid on
xlabel('ka/ka_0','fontsize',20, 'fontangle','italic');
ylabel('TVR [dB]','fontsize',20, 'fontangle','italic');
set(gca, 'fontsize',16)
set(gcf, 'color', 'w')

%%
for i = 1:length(freq)
    Z_11(i,1) = VF_For_FFR{i,1}(1,1);
    Z_12(i,1) = VF_For_FFR{i,1}(1,2);
    Z_13(i,1) = VF_For_FFR{i,1}(1,3);
    Z_22(i,1) = VF_For_FFR{i,1}(2,2);
end
Z_temp = [Z_11 Z_12 Z_13 Z_22];
Z = [real(Z_temp) imag(Z_temp)];
