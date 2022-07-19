TotalLengthofData = length(KaRatio);
gaRatioVector = g_Parameter./a;
gLRatioVector = g_Parameter./L;

x_axis = gaRatioVector;

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

k=1;
for i = 1 : length(MaxTVRLevel)
    FindSignStart = sign(FirstStep(:,i));
    for j = 1 : length(FindSignStart)-1
        if (FindSignStart(j) == -1) && ((FindSignStart(j) + FindSignStart(j+1)) == 0)
            Temp(k) = j+1;
            k=k+1;
        end
    end
    StartPosition(:,i) = max(Temp);
    k=1;
    clear Temp
    
    FindSignEnd = sign(SecondStep(:,i));
    for j = 1 : length(FindSignEnd)-1
        if (FindSignEnd(j) == 1) && ((FindSignEnd(j) + FindSignEnd(j+1)) == 0)
            Temp(k) = j;
            k=k+1;
        end
    end
    EndPositionTemp = min(Temp);
    EndPosition(:,i) = EndPositionTemp + MaxTVRPosition(i) - 1;
    clear Temp 
end

StartKaRatio = KaRatio(StartPosition);
EndKaRatio = KaRatio(EndPosition);

Bandwidth = KaRatio(EndPosition) - KaRatio(StartPosition);

%%
% %% Plot
% figure(1)
% plot(x_axis,StartKaRatio,'LineWidth',2)
% % xlim([0.08 0.2])
% grid on
% title('Start position of Bandwidth','fontsize',24, 'fontangle','italic');
% xlabel('g/L ratio','fontsize',24, 'fontangle','italic');
% ylabel('[ (ka) / (ka)_0 ]','fontsize',24, 'fontangle','italic');
% set(gca, 'fontsize',20)
% set(gcf, 'color', 'w')
% 
% figure(2)
% plot(x_axis,EndKaRatio,'LineWidth',2)
% % xlim([0.08 0.2])
% grid on
% title('End position of Bandwidth','fontsize',24, 'fontangle','italic');
% xlabel('g/L ratio','fontsize',24, 'fontangle','italic');
% ylabel('[ (ka) / (ka)_0 ]','fontsize',24, 'fontangle','italic');
% set(gca, 'fontsize',20)
% set(gcf, 'color', 'w')
% 
% figure(3)
% plot(x_axis,Bandwidth,'LineWidth',2)
% % xlim([0.08 0.2])
% grid on
% title('Bandwidth','fontsize',24, 'fontangle','italic');
% xlabel('g/L ratio','fontsize',24, 'fontangle','italic');
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
% xlabel('g/L ratio','fontsize',24, 'fontangle','italic');
% ylabel('[ (ka) / (ka)_0 ]','fontsize',24, 'fontangle','italic');
% set(gca, 'fontsize',20)
% set(gcf, 'color', 'w')
% 
% subplot(3,1,2)
% plot(x_axis,EndKaRatio,'LineWidth',2)
% % xlim([0.08 0.2])
% grid on
% title('End position of Bandwidth','fontsize',24, 'fontangle','italic');
% xlabel('g/L ratio','fontsize',24, 'fontangle','italic');
% ylabel('[ (ka) / (ka)_0 ]','fontsize',24, 'fontangle','italic');
% set(gca, 'fontsize',20)
% set(gcf, 'color', 'w')
% 
% subplot(3,1,3)
% plot(x_axis,Bandwidth,'LineWidth',2)
% % xlim([0.08 0.2])
% grid on
% title('Bandwidth','fontsize',24, 'fontangle','italic');
% xlabel('g/L ratio','fontsize',24, 'fontangle','italic');
% ylabel('Bandwidth','fontsize',20, 'fontangle','italic');
% set(gca, 'fontsize',20)
% set(gcf, 'color', 'w')
% 
% figure(5)
% plot(x_axis,MaxTVRLevel,'LineWidth',2)
% % xlim([0.08 0.2])
% grid on
% title('TVR','fontsize',24, 'fontangle','italic');
% xlabel('g/L ratio','fontsize',24, 'fontangle','italic');
% ylabel('TVR [dB]','fontsize',24, 'fontangle','italic');
% set(gca, 'fontsize',20)
% set(gcf, 'color', 'w')

%% Plot
figure(1)
plot(x_axis,StartKaRatio,'LineWidth',2)
% xlim([0.08 0.2])
grid on
title('Start position of Bandwidth','fontsize',24, 'fontangle','italic');
xlabel('g/a ratio','fontsize',24, 'fontangle','italic');
ylabel('[ (ka) / (ka)_0 ]','fontsize',24, 'fontangle','italic');
set(gca, 'fontsize',20)
set(gcf, 'color', 'w')

figure(2)
plot(x_axis,EndKaRatio,'LineWidth',2)
% xlim([0.08 0.2])
grid on
title('End position of Bandwidth','fontsize',24, 'fontangle','italic');
xlabel('g/a ratio','fontsize',24, 'fontangle','italic');
ylabel('[ (ka) / (ka)_0 ]','fontsize',24, 'fontangle','italic');
set(gca, 'fontsize',20)
set(gcf, 'color', 'w')

figure(3)
plot(x_axis,Bandwidth,'LineWidth',2)
% xlim([0.08 0.2])
grid on
title('Bandwidth','fontsize',24, 'fontangle','italic');
xlabel('g/a ratio','fontsize',24, 'fontangle','italic');
ylabel('Bandwidth','fontsize',24, 'fontangle','italic');
set(gca, 'fontsize',20)
set(gcf, 'color', 'w')

figure(4)
subplot(3,1,1)
plot(x_axis,StartKaRatio,'LineWidth',2)
% xlim([0.08 0.2])
grid on
title('Start position of Bandwidth','fontsize',24, 'fontangle','italic');
xlabel('g/a ratio','fontsize',24, 'fontangle','italic');
ylabel('[ (ka) / (ka)_0 ]','fontsize',24, 'fontangle','italic');
set(gca, 'fontsize',20)
set(gcf, 'color', 'w')

subplot(3,1,2)
plot(x_axis,EndKaRatio,'LineWidth',2)
% xlim([0.08 0.2])
grid on
title('End position of Bandwidth','fontsize',24, 'fontangle','italic');
xlabel('g/a ratio','fontsize',24, 'fontangle','italic');
ylabel('[ (ka) / (ka)_0 ]','fontsize',24, 'fontangle','italic');
set(gca, 'fontsize',20)
set(gcf, 'color', 'w')

subplot(3,1,3)
plot(x_axis,Bandwidth,'LineWidth',2)
% xlim([0.08 0.2])
grid on
title('Bandwidth','fontsize',24, 'fontangle','italic');
xlabel('g/a ratio','fontsize',24, 'fontangle','italic');
ylabel('Bandwidth','fontsize',20, 'fontangle','italic');
set(gca, 'fontsize',20)
set(gcf, 'color', 'w')

figure(5)
plot(x_axis,MaxTVRLevel,'LineWidth',2)
% xlim([0.08 0.2])
grid on
title('TVR','fontsize',24, 'fontangle','italic');
xlabel('g/a ratio','fontsize',24, 'fontangle','italic');
ylabel('TVR [dB]','fontsize',24, 'fontangle','italic');
set(gca, 'fontsize',20)
set(gcf, 'color', 'w')

%%
% %% Plot
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