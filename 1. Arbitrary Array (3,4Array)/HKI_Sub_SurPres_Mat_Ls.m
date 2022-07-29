%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate HKI for Arbitary Axis-Symmetric Cylinder,
% Result : Surface pressure & RadiationImpedance
% : HKI_Sub_SurPres_Mat_Ls.m
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

function [GD, GW] = HKI_Sub_SurPres_Mat_Ls(freq,RealAcuteAngle,c0,rho0,Line_Sec,NN)

%% Line 정의 및 Normal Vector 검증
totNumOfLine = length(Line_Sec);

Line = zeros(totNumOfLine+1,2);
Line(1:2,:) = Line_Sec(:,:,1);
for LineNum = 2:totNumOfLine
    Line(LineNum+1,:) = Line_Sec(2,:,LineNum);
end

% 각 Line 의 방향 vector 판단
for CalDVec_num = 1:totNumOfLine
    DVector(CalDVec_num,:) = Line(CalDVec_num+1,:) - Line(CalDVec_num,:);
    MagDVector(CalDVec_num,1) = norm(DVector(CalDVec_num,:));
    UnitDVector(CalDVec_num,:) = DVector(CalDVec_num,:)./MagDVector(CalDVec_num);
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
kLmin = (2*pi*freq/c0)*Lmin; % ka 대신 kLim 사용됨

%% 이 부분 개선 필요
TotNumPre = sum(NN)-(totNumOfLine-1); % Total number of Pressure Vector
TotNumVel = sum(NN); % Total number of Velocity Vector

for num = 1:totNumOfLine
    % 각 라인에 노드 배정, IntError 결정, 2%
    detN(num) = nonDimLine_Mag(num) / (NN(num)-1);
    Node_Ls{:,num} = (0:detN(num):nonDimLine_Mag(num))';
    Node_Ls_Ori{:,num} = Node_Ls{:,num};
    IntError(num) = detN(num)*0.002;
    
    % 각 라인의 각도 및 기준 좌표 찾기 및 기준 좌표를 이용하여 전체 형상 노드 변수 생성
    SortLine = sortrows(Line_Sec(:,:,num));
    psi(num) = asin((SortLine(2,2) - SortLine(1,2))/MagDVector(num));
    
    if psi(num) > 0 || psi(num) == 0
        if mean(SortLine(:,2) > 0) == 1
            if (Line_Sec(2,2,num) - Line_Sec(1,2,num)) < 0 || (Line_Sec(2,2,num) - Line_Sec(1,2,num)) == 0
                Node_Ls{:,num} = flipud(Node_Ls{:,num});
                Ref_rzn(num,:) = [SortLine(1,1) SortLine(1,2) SortLine(2,1)];
                Ref_etaztan(num,:)= Ref_rzn(num,:)./Lmin;
                Node_temp(:,1) = Ref_etaztan(num,1) + Node_Ls{:,num}*cos(psi(num));
                Node_temp(:,2) = Ref_etaztan(num,2) + Node_Ls{:,num}*sin(psi(num));
            elseif (Line_Sec(2,2,num) - Line_Sec(1,2,num)) > 0
                Ref_rzn(num,:) = [SortLine(1,1) SortLine(1,2) SortLine(2,1)];
                Ref_etaztan(num,:)= Ref_rzn(num,:)./Lmin;
                Node_temp(:,1) = Ref_etaztan(num,1) + Node_Ls{:,num}*cos(psi(num));
                Node_temp(:,2) = Ref_etaztan(num,2) + Node_Ls{:,num}*sin(psi(num));
            end
        else
            Ref_rzn(num,:) = [SortLine(1,1) SortLine(1,2) SortLine(2,1)];
            Ref_etaztan(num,:)= Ref_rzn(num,:)./Lmin;
            Node_temp(:,1) = Ref_etaztan(num,1) + Node_Ls{:,num}*cos(psi(num));
            Node_temp(:,2) = Ref_etaztan(num,2) + Node_Ls{:,num}*sin(psi(num));
        end
    elseif psi(num) < 0
        if mean(SortLine(:,2) > 0) == 1
            Node_Ls{:,num} = flipud(Node_Ls{:,num});
            Ref_rzn(num,:) = [SortLine(1,1) SortLine(1,2) SortLine(2,1)];
            Ref_etaztan(num,:)= Ref_rzn(num,:)./Lmin;
            Node_temp(:,1) = Ref_etaztan(num,1) + Node_Ls{:,num}*cos(psi(num));
            Node_temp(:,2) = Ref_etaztan(num,2) + Node_Ls{:,num}*sin(psi(num));
        else
            Ref_rzn(num,:) = [SortLine(1,1) SortLine(1,2) SortLine(2,1)];
            Ref_etaztan(num,:)= Ref_rzn(num,:)./Lmin;
            Node_temp(:,1) = Ref_etaztan(num,1) + Node_Ls{:,num}*cos(psi(num));
            Node_temp(:,2) = Ref_etaztan(num,2) + Node_Ls{:,num}*sin(psi(num));
        end
    end
    
    if num == 1
        Node(:,:) = Node_temp;
    else
        Node(sum(NN(1:num-1))+1:sum(NN(1:num)),:) = Node_temp;
    end
    
    clear Node_temp tblA SortTblA
end

for num = 1:(totNumOfLine-1)
    Node(sum(NN(1:num))-(num-1),:) = [];
end

LsNode = Node_Ls; % Cell형 변수, LsNode{ } 형태로 사용
etaNode = Node(:,1);
ztaNode = Node(:,2);

NormalizedDist(1:NN(1),1) = Node_Ls_Ori{1,1}(:,1);
for DistNum = 2:length(Node_Ls_Ori)
    NormalizedDist(length(NormalizedDist)+1:length(NormalizedDist)+NN(DistNum)-1,1) = ...
        NormalizedDist(end,1) + Node_Ls_Ori{1,DistNum}(2:end,1);
end
DimensionDist = NormalizedDist.*Lmin;

%% Solid Angle 정의
for SANum = 1:totNumOfLine+1
    if SANum == 1 %시작 지점
        V1 = NorDVector(1,1:2);
        V2 = [1 0];
        AngleOfNorDV(1,1) = round(asind((V1(1)*V2(2)-V1(2)*V2(1))/(norm(V1)*norm(V2))));
    elseif SANum == totNumOfLine+1 % 끝지점
        V2 = NorDVector(totNumOfLine,1:2);
        V1 = [1 0];
        AngleOfNorDV(totNumOfLine+1,1) = round(asind((V1(1)*V2(2)-V1(2)*V2(1))/(norm(V1)*norm(V2))));
    else % 나머지 부분
        V1 = NorDVector(SANum-1,1:2);
        V2 = NorDVector(SANum,1:2);
        AngleOfNorDV(SANum,1) = round(asind((V1(1)*V2(2)-V1(2)*V2(1))/(norm(V1)*norm(V2))));
    end

    if RealAcuteAngle(SANum,1) == 0
        if SANum == 1 || SANum == totNumOfLine+1 %시작 지점, 끝지점
            Angle(SANum,1) = (180-(AngleOfNorDV(SANum,1)));
        else
            Angle(SANum,1) = (180+AngleOfNorDV(SANum,1));
        end
    elseif RealAcuteAngle(SANum,1) == 1
        if SANum == 1 || SANum == totNumOfLine+1 %시작 지점, 끝지점
            Angle(SANum,1) = AngleOfNorDV(SANum,1);
        else
            if AngleOfNorDV(SANum,1) < 0
                Angle(SANum,1) = abs(AngleOfNorDV(SANum,1));
            else
                Angle(SANum,1) = (360-AngleOfNorDV(SANum,1));
            end
        end
    end
    
    if SANum == 1 || SANum == totNumOfLine+1  %시작 지점, 끝지점
        SAatEdge(SANum,1) = 2*pi.*(1-cosd(Angle(SANum,1)));
    else
        SAatEdge(SANum,1) = (Angle(SANum,1)/360)*4*pi; % 각 Line이 만나는 부분의 각도
    end
end
clear SANum

SolidAngle = ones(length(etaNode),1)*2*pi;
SolidAngle(1,1) = SAatEdge(1,1);
SolidAngle(end,1) = SAatEdge(end,1);

AlReadySum = 1;
for num = 2:totNumOfLine
    SANum = AlReadySum + (NN(num-1)-1);
    SolidAngle(SANum,1) = SAatEdge(num,1);
    AlReadySum = SANum;
end

%% Cal Global D Matrix
% GD(PressureVectorNumber, VelocityVectorNumber)
GD = zeros(TotNumPre,TotNumVel);
MD = cell(1,totNumOfLine);
ThetaSliceNum = 9;
IntTheta = linspace(0,2*pi,ThetaSliceNum);

for LineNum = 1 : totNumOfLine
    CaseNum = NorDVector(LineNum,3);
    
    Loacl_NorDVec = NorDVector(LineNum,1:2);
    Local_psi = psi(LineNum);
    Local_Ref_etaztan = Ref_etaztan(LineNum,:);
    Local_LineN = NN(LineNum);
    Local_LsNode = LsNode{:,LineNum};
    Local_detN = detN(LineNum);
    Local_IntError = IntError(LineNum);
    
    parfor num = 1:TotNumPre
        % Define eto_o, zta_o
        eta_o = etaNode(num);
        zta_o = ztaNode(num);
        
        % Define Solid Angle
        Local_SolidAngle = SolidAngle(num);

        R=@(eta,theta,zta) HKI_Sub_R(eta_o,zta_o,eta,theta,zta); % Define R using Sub_R(eta_o,zta_o,eta,theta,zta,LaRatio)
        Fac=((1i*kLmin)/Local_SolidAngle); % After Integral, multiply Fac.
        
        MD_1 = zeros(1,Local_LineN); % Temp Variable; MultiThread
        MD_2 = zeros(1,Local_LineN); % Temp Variable; MultiThread
        
        for num_b = 1:Local_LineN
            IntRange = [Local_LsNode(num_b)-Local_detN Local_LsNode(num_b)];
            Int_D=@(Ls,theta) HKI_Sub_IntFn_D_Ls(CaseNum,R,kLmin,Local_psi,Local_Ref_etaztan,IntRange,Loacl_NorDVec,Ls,theta); % Define Integral Equation using Sub_IntFn_D(Surface Num,R,ka,eta,theta,zta)
            Int_Fn=@(Ls,theta) HKI_Sub_SFn(1,Local_LsNode,Ls,IntRange).*Int_D(Ls,theta); % Define Integral Fn, Sub_SFn is define Shape Function
            
%             for num_th = 1:(length(IntTheta)-1)
%                 MD_1(num_b) = MD_1(num_b) + HKI_Sub_Integral_Ls(Int_Fn,IntRange,IntTheta(num_th),IntTheta(num_th+1),Local_IntError)*Fac;
%             end
                
            MD_1(num_b) = HKI_Sub_Integral_Ls_theta(Int_Fn,IntRange,IntTheta,Local_IntError)*Fac;
            
            IntRange =  [Local_LsNode(num_b) Local_LsNode(num_b)+Local_detN];
            Int_D=@(Ls,theta) HKI_Sub_IntFn_D_Ls(CaseNum,R,kLmin,Local_psi,Local_Ref_etaztan,IntRange,Loacl_NorDVec,Ls,theta); % Define Integral Equation using Sub_IntFn_D(Surface Num,R,ka,eta,theta,zta)
            Int_Fn=@(Ls,theta) HKI_Sub_SFn(2,Local_LsNode,Ls,IntRange).*Int_D(Ls,theta); % Define Integral Fn
            
%             for num_th = 1:(length(IntTheta)-1)
%                 MD_2(num_b) = MD_2(num_b) + HKI_Sub_Integral_Ls(Int_Fn,IntRange,IntTheta(num_th),IntTheta(num_th+1),Local_IntError)*Fac;
%             end
            
            MD_2(num_b) = HKI_Sub_Integral_Ls_theta(Int_Fn,IntRange,IntTheta,Local_IntError)*Fac;
        end
        MD_3(num,:) = MD_1 + MD_2;
    end
    MD{1,LineNum} = MD_3; % MD matrix construction; MultiThread
    clear MD_3
end
clear Loacl_NorDVec Local_psi Local_Ref_etaztan Local_LineN Local_LsNode Local_detN Local_IntError

% GD = zeros(TotNumPre,TotNumVel);
% MD = cell(1,totNumOfLine);
% 
% for LineNum = 1 : totNumOfLine
%     CaseNum = NorDVector(LineNum,3);
%     
%     Loacl_NorDVec = NorDVector(LineNum,1:2);
%     Local_psi = psi(LineNum);
%     Local_Ref_etaztan = Ref_etaztan(LineNum,:);
%     Local_LineN = NN(LineNum);
%     Local_LsNode = LsNode{:,LineNum};
%     Local_detN = detN(LineNum);
%     Local_IntError = IntError(LineNum);
%     
%     parfor num = 1:TotNumPre
%         % Define eto_o, zta_o
%         eta_o = etaNode(num);
%         zta_o = ztaNode(num);
%         
%         % Define Solid Angle
%         Local_SolidAngle = SolidAngle(num);
% 
%         R=@(eta,theta,zta) HKI_Sub_R(eta_o,zta_o,eta,theta,zta); % Define R using Sub_R(eta_o,zta_o,eta,theta,zta,LaRatio)
%         Fac=((1i*kLmin)/Local_SolidAngle); % After Integral, multiply Fac.
%         
%         MD_1 = zeros(1,Local_LineN); % Temp Variable; MultiThread
%         MD_2 = zeros(1,Local_LineN); % Temp Variable; MultiThread
%         
%         for num_b = 1:Local_LineN
%             IntRange = [Local_LsNode(num_b)-Local_detN Local_LsNode(num_b)];
%             Int_D=@(Ls,theta) HKI_Sub_IntFn_D_Ls(CaseNum,R,kLmin,Local_psi,Local_Ref_etaztan,IntRange,Loacl_NorDVec,Ls,theta); % Define Integral Equation using Sub_IntFn_D(Surface Num,R,ka,eta,theta,zta)
%             Int_Fn=@(Ls,theta) HKI_Sub_SFn(1,Local_LsNode,Ls,IntRange).*Int_D(Ls,theta); % Define Integral Fn, Sub_SFn is define Shape Function
%             MD_1(num_b) = HKI_Sub_Integral(Int_Fn,IntRange,Local_IntError)*Fac;
%             
%             IntRange =  [Local_LsNode(num_b) Local_LsNode(num_b)+Local_detN];
%             Int_D=@(Ls,theta) HKI_Sub_IntFn_D_Ls(CaseNum,R,kLmin,Local_psi,Local_Ref_etaztan,IntRange,Loacl_NorDVec,Ls,theta); % Define Integral Equation using Sub_IntFn_D(Surface Num,R,ka,eta,theta,zta)
%             Int_Fn=@(Ls,theta) HKI_Sub_SFn(2,Local_LsNode,Ls,IntRange).*Int_D(Ls,theta); % Define Integral Fn
%             MD_2(num_b) = HKI_Sub_Integral(Int_Fn,IntRange,Local_IntError)*Fac;
%         end
%         MD_3(num,:) = MD_1 + MD_2; % MD matrix construction; MultiThread
%     end
%     MD{1,LineNum} = MD_3; % MD matrix construction; MultiThread
%     clear MD_3
% end
% clear Loacl_NorDVec Local_psi Local_Ref_etaztan Local_LineN Local_LsNode Local_detN Local_IntError

%% Cal Global W Matrix
GW = zeros(TotNumPre,TotNumPre);
GW_UnSum = zeros(TotNumPre,TotNumVel); % UnSum
MW = cell(1,totNumOfLine);

for LineNum = 1 : totNumOfLine
    CaseNum = NorDVector(LineNum,3);

    Loacl_NorDVec = NorDVector(LineNum,1:2);
    Local_psi = psi(LineNum);
    Local_Ref_etaztan = Ref_etaztan(LineNum,:);
    Local_LineN = NN(LineNum);
    Local_LsNode = LsNode{:,LineNum};
    Local_detN = detN(LineNum);
    Local_IntError = IntError(LineNum);
    
    parfor num = 1:TotNumPre
        % Define eto_o, zta_o
        eta_o = etaNode(num);
        zta_o = ztaNode(num);
        
        % Define Solid Angle
        Local_SolidAngle = SolidAngle(num);
        
        R=@(eta,theta,zta) HKI_Sub_R(eta_o,zta_o,eta,theta,zta); % Define R
        Fac=1/(Local_SolidAngle);
        
        MW_1 = zeros(1,Local_LineN); % Temp Variable; MultiThread
        MW_2 = zeros(1,Local_LineN); % Temp Variable; MultiThread
        
        for num_b = 1:Local_LineN
            IntRange = [Local_LsNode(num_b)-Local_detN Local_LsNode(num_b)];
            Int_W=@(Ls,theta) HKI_Sub_IntFn_W_Ls(CaseNum,R,kLmin,Local_psi,Local_Ref_etaztan,eta_o,zta_o,IntRange,Loacl_NorDVec,Ls,theta); % Define Integral Equation
            Int_Fn=@(Ls,theta) HKI_Sub_SFn(1,Local_LsNode,Ls,IntRange).*Int_W(Ls,theta); % Define Integral Fn
            
%             for num_th = 1:(length(IntTheta)-1)
%                 MW_1(num_b) = MW_1(num_b) + HKI_Sub_Integral_Ls(Int_Fn,IntRange,IntTheta(num_th),IntTheta(num_th+1),Local_IntError)*Fac;
%             end
            
            MW_1(num_b) = HKI_Sub_Integral_Ls_theta(Int_Fn,IntRange,IntTheta,Local_IntError)*Fac;
            
            IntRange =  [Local_LsNode(num_b) Local_LsNode(num_b)+Local_detN];
            Int_W=@(Ls,theta) HKI_Sub_IntFn_W_Ls(CaseNum,R,kLmin,Local_psi,Local_Ref_etaztan,eta_o,zta_o,IntRange,Loacl_NorDVec,Ls,theta); % Define Integral Equation
            Int_Fn=@(Ls,theta) HKI_Sub_SFn(2,Local_LsNode,Ls,IntRange).*Int_W(Ls,theta); % Define Integral Fn
            
%             for num_th = 1:(length(IntTheta)-1)
%                 MW_2(num_b) = MW_2(num_b) + HKI_Sub_Integral_Ls(Int_Fn,IntRange,IntTheta(num_th),IntTheta(num_th+1),Local_IntError)*Fac;
%             end
            
            MW_2(num_b) = HKI_Sub_Integral_Ls_theta(Int_Fn,IntRange,IntTheta,Local_IntError)*Fac;
        end
        MW_3(num,:) = MW_1 + MW_2; % MW matrix construction; MultiThread
    end
    MW{1,LineNum} = MW_3; % MD matrix construction; MultiThread
    clear MW_3
end

% GW = zeros(TotNumPre,TotNumPre);
% GW_UnSum = zeros(TotNumPre,TotNumVel); % UnSum
% MW = cell(1,totNumOfLine);
% 
% for LineNum = 1 : totNumOfLine
%     CaseNum = NorDVector(LineNum,3);
% 
%     Loacl_NorDVec = NorDVector(LineNum,1:2);
%     Local_psi = psi(LineNum);
%     Local_Ref_etaztan = Ref_etaztan(LineNum,:);
%     Local_LineN = NN(LineNum);
%     Local_LsNode = LsNode{:,LineNum};
%     Local_detN = detN(LineNum);
%     Local_IntError = IntError(LineNum);
%     
%     parfor num = 1:TotNumPre
%         % Define eto_o, zta_o
%         eta_o = etaNode(num);
%         zta_o = ztaNode(num);
%         
%         % Define Solid Angle
%         Local_SolidAngle = SolidAngle(num);
%         
%         R=@(eta,theta,zta) HKI_Sub_R(eta_o,zta_o,eta,theta,zta); % Define R
%         Fac=1/(Local_SolidAngle);
%         
%         MW_1 = zeros(1,Local_LineN); % Temp Variable; MultiThread
%         MW_2 = zeros(1,Local_LineN); % Temp Variable; MultiThread
%         
%         for num_b = 1:Local_LineN
%             IntRange = [Local_LsNode(num_b)-Local_detN Local_LsNode(num_b)];
%             Int_W=@(Ls,theta) HKI_Sub_IntFn_W_Ls(CaseNum,R,kLmin,Local_psi,Local_Ref_etaztan,eta_o,zta_o,IntRange,Loacl_NorDVec,Ls,theta); % Define Integral Equation
%             Int_Fn=@(Ls,theta) HKI_Sub_SFn(1,Local_LsNode,Ls,IntRange).*Int_W(Ls,theta); % Define Integral Fn
%                        
%             MW_1(num_b) = HKI_Sub_Integral(Int_Fn,IntRange,Local_IntError)*Fac;
%             
%             IntRange =  [Local_LsNode(num_b) Local_LsNode(num_b)+Local_detN];
%             Int_W=@(Ls,theta) HKI_Sub_IntFn_W_Ls(CaseNum,R,kLmin,Local_psi,Local_Ref_etaztan,eta_o,zta_o,IntRange,Loacl_NorDVec,Ls,theta); % Define Integral Equation
%             Int_Fn=@(Ls,theta) HKI_Sub_SFn(2,Local_LsNode,Ls,IntRange).*Int_W(Ls,theta); % Define Integral Fn
% 
%             MW_2(num_b) = HKI_Sub_Integral(Int_Fn,IntRange,Local_IntError)*Fac;
%         end
%         MW_3(num,:) = MW_1 + MW_2; % MW matrix construction; MultiThread
%     end
%     MW{1,LineNum} = MW_3; % MD matrix construction; MultiThread
%     clear MW_3
% end

%% Matrix Calculation
AlReadySum = 1;
for LineNum = 1:totNumOfLine
    GDNumStart = AlReadySum;
    GDNumEnd = (AlReadySum + (NN(LineNum)))-1;
    GD(:,GDNumStart:GDNumEnd) = MD{1,LineNum};
    GW_UnSum(:,GDNumStart:GDNumEnd) = MW{1,LineNum};
    AlReadySum = GDNumEnd+1;
end

AlReadySum = 1;
GW(:,1) = MW{1,1}(:,1);
for LineNum = 1:totNumOfLine
    GWNumStart = AlReadySum+1;
    GWNumEnd = (AlReadySum + (NN(LineNum)))-2;
    GW(:,GWNumStart:GWNumEnd) = MW{1,LineNum}(:,2:NN(LineNum)-1);
    if LineNum ~= totNumOfLine
        GW(:,GWNumEnd+1) = MW{1,LineNum}(:,NN(LineNum))+MW{1,LineNum+1}(:,1);
    elseif LineNum == totNumOfLine
        GW(:,GWNumEnd+1) = MW{1,LineNum}(:,NN(LineNum));
    end
    AlReadySum = GWNumEnd+1;
end

end