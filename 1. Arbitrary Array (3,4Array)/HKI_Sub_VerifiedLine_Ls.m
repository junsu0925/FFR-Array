function [VerifiedLine,RealAcuteAngle,MagDVector] = HKI_Sub_VerifiedLine_Ls(Line_Sec)
totNumOfLine = length(Line_Sec);
RealAcuteAngle = zeros(totNumOfLine+1,1);

disp(' ')
disp('First, we will verify the information of each input line and the calculated vertical direction vector result.')

Line = zeros(totNumOfLine+1,2);
Line(1:2,:) = Line_Sec(:,:,1);
for LineNum = 2:totNumOfLine
    Line(LineNum+1,:) = Line_Sec(2,:,LineNum);
end

% 각 Line 의 방향 vector 판단
for CalDVec_num = 1:totNumOfLine
    DVector = Line(CalDVec_num+1,:) - Line(CalDVec_num,:);
    MagDVector(CalDVec_num) = norm(DVector);
    UnitDVector(CalDVec_num,:)  = DVector./MagDVector(CalDVec_num);
    NorDVector(CalDVec_num,1:2) = [UnitDVector(CalDVec_num,2) -UnitDVector(CalDVec_num,1)];
    
    if abs(NorDVector(CalDVec_num,1)) == 1
        NorDVector(CalDVec_num,3) = 2;
    elseif abs(NorDVector(CalDVec_num,2)) == 1
        NorDVector(CalDVec_num,3) = 1;
    else
        NorDVector(CalDVec_num,3) = 3;
    end    
end
clear CalDVec_num
Lmin = min(MagDVector); % Dimensionless 를 위한 기준 길이 선정
nonDimLine_Mag = MagDVector/Lmin; % 각 Line 길이 Dimensionless 계산

%% Vector plot 출력, 확인
LineCen = zeros(totNumOfLine,2);
for LineCenNum = 1:totNumOfLine
    if NorDVector(LineCenNum,3) == 1
        LineCen(LineCenNum,:) = [(Line_Sec(1,1,LineCenNum)+Line_Sec(2,1,LineCenNum))/2 Line_Sec(1,2,LineCenNum)];
    elseif NorDVector(LineCenNum,3) == 2
        LineCen(LineCenNum,:) = [Line_Sec(1,1,LineCenNum) (Line_Sec(1,2,LineCenNum)+Line_Sec(2,2,LineCenNum))/2];
    elseif NorDVector(LineCenNum,3) == 3
        LineCen(LineCenNum,:) = [(Line_Sec(1,1,LineCenNum)+Line_Sec(2,1,LineCenNum))/2 ...
            (Line_Sec(1,2,LineCenNum)+Line_Sec(2,2,LineCenNum))/2];
    end
end
clear LineCenNum
fig1 = figure(9999);
rmax = max(Line(:,1));
zmax = (max(Line(:,2))-min(Line(:,2)));
set(fig1,'units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
plot(Line(:,1),Line(:,2),'LineWidth', 2)
hold on
quiver(LineCen(:,1),LineCen(:,2),NorDVector(:,1)*Lmin/3,NorDVector(:,2)*Lmin/3,'AutoScale','off','LineWidth', 2)
grid on
pbaspect([rmax zmax 1])
xlabel('r [m]','fontsize',20, 'fontangle','italic');
ylabel('z [m]','fontsize',20, 'fontangle','italic');
legend('Geometry','Outer Normal Vector','Location', 'Best')
set(gca, 'fontsize',16)
set(gcf, 'color', 'w')

pause(0.1)
disp(' ')
Conti = input('Are the geometry graphs drawn as desired? (Yes = 1, No = 2) : ');
if Conti == 1
    VerifiedLine = 1;
    %% Angle Check
    for SANum = 1:totNumOfLine+1
        if SANum == 1 %시작 지점
            V2 = NorDVector(1,1:2);
            V1 = [-V2(1) V2(2)];
            AngleOfNorDV(1,1) = round(asind((V1(1)*V2(2)-V1(2)*V2(1))/(norm(V1)*norm(V2))));
        elseif SANum == totNumOfLine+1 % 끝지점
            V1 = NorDVector(totNumOfLine,1:2);
            V2 = [-V1(1) V1(2)];
            AngleOfNorDV(totNumOfLine+1,1) = round(asind((V1(1)*V2(2)-V1(2)*V2(1))/(norm(V1)*norm(V2))));
        else % 나머지 부분
            V1 = NorDVector(SANum-1,1:2);
            V2 = NorDVector(SANum,1:2);
            AngleOfNorDV(SANum,1) = round(asind((V1(1)*V2(2)-V1(2)*V2(1))/(norm(V1)*norm(V2))));
        end
    end
    
    FindAcuteAngle = find(abs(AngleOfNorDV)<90 & AngleOfNorDV ~=0);
    if isempty(FindAcuteAngle) == 0
        RealAcuteAngle = zeros(length(AngleOfNorDV),1);
        disp(' ')
        disp('As a result of calculating the angle with the input coordinates,')
        disp('it seems that acute angle exists. Which point is really acute?')
        disp(' ')
        disp('The expected location of the acute angle is as follows.')
        AcuteAngle = FindAcuteAngle;
        tblA = table(AcuteAngle);
        disp(' ')
        disp(tblA)
        disp(' ')
        RealAcuteAngle_temp = input('Real Acute Angle position, If not, just Enter ex, [1 2 3] : ');
        RealAcuteAngle(RealAcuteAngle_temp',1) = 1;
    end
elseif Conti == 2
    VerifiedLine = 2;
    disp(' ')
    disp('This code will be terminated soon.')
end
end