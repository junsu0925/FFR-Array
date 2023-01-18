KaRatio = Wave.KaRatio(:,1);
TotalLengthofData = length(KaRatio);

gaRatioVector = Geometry.gSideValue./Geometry.a;
gLRatioVector = Geometry.gSideValue./Geometry.RL(1);

x_axis = gLRatioVector;

[MaxTVRLevel,MaxTVRPosition] = max(TVRSave);
MinTVRLevel = MaxTVRLevel - 6;

%% Find Bandwidth
FirstStep = ones(max(MaxTVRPosition),length(MaxTVRLevel));
for i = 1 : length(MaxTVRLevel)
    Temp = TVRSave(1:MaxTVRPosition(i),i) - MinTVRLevel(i);
    FirstStep(1:length(Temp),i) = Temp;
end
clear Temp

SecondStep = ones((TotalLengthofData-min(MaxTVRPosition))+1,length(MaxTVRLevel));
for i = 1 : length(MaxTVRLevel)
    Temp = TVRSave(MaxTVRPosition(i):end,i) - MinTVRLevel(i);
    SecondStep(1:length(Temp),i) = Temp;
end
clear Temp

for i = 1 : length(MaxTVRLevel)
    [StartDiffLevel(:,i),StartPosition(:,i)]=min(abs(FirstStep(:,i)));
    [EndDiffLevel(:,i),EndPositionTemp]=min(abs(SecondStep(:,i)));
    EndPosition(:,i) = EndPositionTemp + MaxTVRPosition(i) - 1;
end

StartKaRatio = KaRatio(StartPosition);
EndKaRatio = KaRatio(EndPosition);

Bandwidth = KaRatio(EndPosition) - KaRatio(StartPosition);

%%
figure(4)
subplot(3,1,1)
plot(x_axis,StartKaRatio,'LineWidth',2)
% xlim([0.08 0.2])
grid on
title('Start position of Bandwidth','fontsize',24, 'fontangle','italic');
xlabel('g/L ratio','fontsize',24, 'fontangle','italic');
ylabel('[ (ka) / (ka)_0 ]','fontsize',24, 'fontangle','italic');
set(gca, 'fontsize',20)
set(gcf, 'color', 'w')

subplot(3,1,2)
plot(x_axis,EndKaRatio,'LineWidth',2)
% xlim([0.08 0.2])
grid on
title('End position of Bandwidth','fontsize',24, 'fontangle','italic');
xlabel('g/L ratio','fontsize',24, 'fontangle','italic');
ylabel('[ (ka) / (ka)_0 ]','fontsize',24, 'fontangle','italic');
set(gca, 'fontsize',20)
set(gcf, 'color', 'w')

subplot(3,1,3)
plot(x_axis,Bandwidth,'LineWidth',2)
% xlim([0.08 0.2])
grid on
title('Bandwidth','fontsize',24, 'fontangle','italic');
xlabel('g/L ratio','fontsize',24, 'fontangle','italic');
ylabel('Bandwidth','fontsize',20, 'fontangle','italic');
set(gca, 'fontsize',20)
set(gcf, 'color', 'w')

figure(5)
plot(x_axis,MaxTVRLevel,'LineWidth',2)
% xlim([0.08 0.2])
grid on
title('TVR','fontsize',24, 'fontangle','italic');
xlabel('g/L ratio','fontsize',24, 'fontangle','italic');
ylabel('TVR [dB]','fontsize',24, 'fontangle','italic');
set(gca, 'fontsize',20)
set(gcf, 'color', 'w')

%%
% figure(4)
% subplot(3,1,1)
% plot(x_axis,StartKaRatio,'LineWidth',2)
% % xlim([0.08 0.2])
% grid on
% title('Start position of Bandwidth','fontsize',24, 'fontangle','italic');
% xlabel('g/a ratio','fontsize',24, 'fontangle','italic');
% ylabel('[ (ka) / (ka)_0 ]','fontsize',24, 'fontangle','italic');
% set(gca, 'fontsize',20)
% set(gcf, 'color', 'w')
% 
% subplot(3,1,2)
% plot(x_axis,EndKaRatio,'LineWidth',2)
% % xlim([0.08 0.2])
% grid on
% title('End position of Bandwidth','fontsize',24, 'fontangle','italic');
% xlabel('g/a ratio','fontsize',24, 'fontangle','italic');
% ylabel('[ (ka) / (ka)_0 ]','fontsize',24, 'fontangle','italic');
% set(gca, 'fontsize',20)
% set(gcf, 'color', 'w')
% 
% subplot(3,1,3)
% plot(x_axis,Bandwidth,'LineWidth',2)
% % xlim([0.08 0.2])
% grid on
% title('Bandwidth','fontsize',24, 'fontangle','italic');
% xlabel('g/a ratio','fontsize',24, 'fontangle','italic');
% ylabel('Bandwidth','fontsize',20, 'fontangle','italic');
% set(gca, 'fontsize',20)
% set(gcf, 'color', 'w')
% 
% figure(5)
% plot(x_axis,MaxTVRLevel,'LineWidth',2)
% % xlim([0.08 0.2])
% grid on
% title('TVR','fontsize',24, 'fontangle','italic');
% xlabel('g/a ratio','fontsize',24, 'fontangle','italic');
% ylabel('TVR [dB]','fontsize',24, 'fontangle','italic');
% set(gca, 'fontsize',20)
% set(gcf, 'color', 'w')

%% Plot
% figure(1)
% plot(x_axis,StartKaRatio,'LineWidth',2)
% % xlim([0.08 0.2])
% grid on
% title('Start position of Bandwidth','fontsize',24, 'fontangle','italic');
% xlabel('Center gap value [m]','fontsize',24, 'fontangle','italic');
% ylabel('[ (ka) / (ka)_0 ]','fontsize',24, 'fontangle','italic');
% set(gca, 'fontsize',20)
% set(gcf, 'color', 'w')
% 
% figure(2)
% plot(x_axis,EndKaRatio,'LineWidth',2)
% % xlim([0.08 0.2])
% grid on
% title('End position of Bandwidth','fontsize',24, 'fontangle','italic');
% xlabel('Center gap value [m]','fontsize',24, 'fontangle','italic');
% ylabel('[ (ka) / (ka)_0 ]','fontsize',24, 'fontangle','italic');
% set(gca, 'fontsize',20)
% set(gcf, 'color', 'w')
% 
% figure(3)
% plot(x_axis,Bandwidth,'LineWidth',2)
% % xlim([0.08 0.2])
% grid on
% title('Bandwidth','fontsize',24, 'fontangle','italic');
% xlabel('Center gap value [m]','fontsize',24, 'fontangle','italic');
% ylabel('Bandwidth','fontsize',24, 'fontangle','italic');
% set(gca, 'fontsize',20)
% set(gcf, 'color', 'w')
% 
% figure(4)
% subplot(3,1,1)
% plot(x_axis,StartKaRatio,'LineWidth',2)
% % xlim([0.08 0.2])
% grid on
% title('Start position of Bandwidth','fontsize',24, 'fontangle','italic');
% xlabel('Center gap value [m]','fontsize',24, 'fontangle','italic');
% ylabel('[ (ka) / (ka)_0 ]','fontsize',24, 'fontangle','italic');
% set(gca, 'fontsize',20)
% set(gcf, 'color', 'w')
% 
% subplot(3,1,2)
% plot(x_axis,EndKaRatio,'LineWidth',2)
% % xlim([0.08 0.2])
% grid on
% title('End position of Bandwidth','fontsize',24, 'fontangle','italic');
% xlabel('Center gap value [m]','fontsize',20, 'fontangle','italic');
% ylabel('[ (ka) / (ka)_0 ]','fontsize',24, 'fontangle','italic');
% set(gca, 'fontsize',20)
% set(gcf, 'color', 'w')
% 
% subplot(3,1,3)
% plot(x_axis,Bandwidth,'LineWidth',2)
% % xlim([0.08 0.2])
% grid on
% title('Bandwidth','fontsize',24, 'fontangle','italic');
% xlabel('Center gap value [m]','fontsize',24, 'fontangle','italic');
% ylabel('Bandwidth','fontsize',20, 'fontangle','italic');
% set(gca, 'fontsize',20)
% set(gcf, 'color', 'w')
% 
% figure(5)
% plot(x_axis,MaxTVRLevel,'LineWidth',2)
% % xlim([0.08 0.2])
% grid on
% title('TVR','fontsize',24, 'fontangle','italic');
% xlabel('Center gap value [m]','fontsize',24, 'fontangle','italic');
% ylabel('TVR [dB]','fontsize',24, 'fontangle','italic');
% set(gca, 'fontsize',20)
% set(gcf, 'color', 'w')