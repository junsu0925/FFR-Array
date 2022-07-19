% Find freq = 1680 : 2410

[n1,m1]=find(f == 1680);
[n2,m2]=find(f == 2410);

TotalNum = (m2-m1)+1;
RealY1Save = real(Y1Save);

for GapNum = 1:gpnum
    for FindNum = m1:m2-1
        MinFirst = RealY1Save(FindNum,GapNum);
        MinSecond = RealY1Save(FindNum+1,GapNum);
        if MinFirst < MinSecond
            ConductanceMin(GapNum,1) = MinFirst;
            break
        end
    end
end

gaRatioVector = g_Parameter./a;
gLRatioVector = g_Parameter./L;

figure(1)
plot(gaRatioVector,ConductanceMin,'LineWidth', 2)
grid on
xlabel('g/a ratio','fontsize',20, 'fontangle','italic');
ylabel('Minimum of Conductance','fontsize',20, 'fontangle','italic');
set(gca, 'fontsize',16)
set(gcf, 'color', 'w')

figure(2)
plot(gLRatioVector,ConductanceMin,'LineWidth', 2)
grid on
xlabel('g/L ratio','fontsize',20, 'fontangle','italic');
ylabel('Minimum of Conductance','fontsize',20, 'fontangle','italic');
set(gca, 'fontsize',16)
set(gcf, 'color', 'w')

    