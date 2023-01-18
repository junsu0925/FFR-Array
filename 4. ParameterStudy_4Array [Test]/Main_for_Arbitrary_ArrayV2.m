%% FFR Circuit Model Code
% Made by Seungwon Nam
% Modified by Kyounghun Been

clc
clear
close all

tic

%% Input Geometry and Material Inputs of FFR
% Assume Matrial of PZT is PZT-5H and Medium is Water

% PZT Transducetion Parameters (PZT-5H)
sE11=1.65e-11;
eT33=3400*8.8541878176*10^-12;
d31=-2.74e-10;
k31=sqrt(d31^2/(sE11*eT33));
eS33=eT33*(1-k31^2);

rhopzt=7500; % Density of PZT
rhowater=999; % Density of Water
sp=1500; % Sound Speed at Water

f=2.8060e+03*(0.3:0.005:1.3);
% f=1000:100:4000; % Set Frequency Range
% f=3000:20:4000; % Set Frequency Range
% f=25:5:4000; % Set Frequency Range
% f=700:5:4000; % Set Frequency Range
% ff=700:5:4000; % Ser Frequency Range for Table

omega=2*pi*f; % Angular Frequency
k=@(w) w/sp; % Wave Number

Geometry.a=0.175; % Ring Average Radius
L=0.1925; % Ring Height
Geometry.t=0.04; % Ring Thickness

% a=0.20; % Ring Average Radius
% L=0.18; % Ring Height
% t=0.04; % Ring Thickness

ka = ((2*pi.*f)./sp).*Geometry.a; % ka, dimensionless parameter by Been
% kaForTable = ((2*pi.*ff)./sp).*a;

% num=input('Number of FFR = '); % Number of FFR
num = 4;
zlength1=4*num-1; % length of the matrix of array number = N
zlength2=2*num+1; % lengths of the matrix of array number = N

R0=2; % the r-distance to calculate TVR
Z0=0; % the z-distance to calculate TVR
% zlength=3+4*(num-1); % length of the impedance matrix of array number = num
zlength=(4*num)-1; % length of the impedance matrix of array number = num

Geometry.RL = zeros(num,1); % Ring Height
for Count = 1 : num
    Geometry.RL(Count,1) = L;
end
clear Count

if num==1 % Single FFR
    Geometry.g=0;
else % Array FFR
    Geometry.g = zeros(num-1,1); % Gab Height
    for Count = 1 : num - 1
        Geometry.g(Count,1) = 0.08;
    end
    clear Count
end

% g_Parameter = 0.08; 
% g_Parameter_Side = [0.02 0.08:0.02:0.2]; 
g_Parameter_Center = [0.02 0.08:0.02:0.2]; 

%% Parameter Study Start

% YSave = zeros(length(f),length(g_Parameter));
% TVRSave = zeros(length(f),length(g_Parameter));
% KaRatioSave = zeros(length(f),length(g_Parameter));

% for gpnum = 1:length(g_Parameter_Side)
for gpnum = 1:length(g_Parameter_Center)
%     Geometry.g(1,1) = g_Parameter_Side(gpnum);
%     Geometry.g(3,1) = Geometry.g(1,1);
    Geometry.g(2,1) = g_Parameter_Center(gpnum);
        
    gaRatio = Geometry.g/Geometry.a; % g/a ratio, dimensionless parameter by Been
    LaRatio = L/Geometry.a; % L/a ratio, dimensionless parameter by Been
    Dimfactor = rhowater*sp*2*pi*Geometry.a*L;
    DimfactorGap = rhowater*sp*2*pi*Geometry.a*Geometry.g(1,1);
    
    Process = sprintf('Parameter g = %1.3f [m]',Geometry.g(1,1))
    
    n=100; % Sum Number (n>50 is enough)    
          
    Geometry.L = sum(Geometry.RL) + sum(Geometry.g); % Total length of FFR array 
    ka = ((2*pi.*f)./sp).*Geometry.a; % ka, dimensionless parameter by Been
    Geometry.LaRatio = Geometry.RL./Geometry.a; % L/a ratio, dimensionless parameter by Been
    Geometry.gaRatio = Geometry.g./Geometry.a; % g/a ratio, dimensionless parameter by Been
    
    n=500; % Sum Number (n>50 is enough)
    
    % Define node m (location)
    Geometry.Zam = zeros(zlength2 - 1,1);
    for Count1 = 1 : zlength2 - 1
        Lai = 0;
        for Count2 = 1 : Count1 - 1
            if mod(Count2,2) == 1 % Odd
                Lai_temp = Geometry.RL((Count2+1)/2);
            elseif mod(Count2,2) == 0 % Even
                Lai_temp = Geometry.g(Count2/2);
            end
            Lai = Lai + Lai_temp;
        end
        Geometry.Zam(Count1,1) = Geometry.L/2 - Lai;
    end
    clear Lai Lai_temp Count1 Count2
    % Define node m (location)
    
    [Geometry.Line_Sec,Geometry.RealAcuteAngle,Geometry.MagDVector] = HKI_Sub_Geometry(Geometry,num);


    %% Circuit Parameters of PZT
    % Ignore Electrical Dissipation G0 and Mechanical Damping Rm
    
    C0=(2*pi*Geometry.a*L/Geometry.t)*eS33; % Equivalent Clamped Capacitance
    N=(2*pi*L)*(d31/sE11); % Equivalent Electromechanical turns Ratio
    M=rhopzt*2*pi*Geometry.a*L*Geometry.t; % Mass of PZT
    CE=(sE11*Geometry.a)/(2*pi*Geometry.t*L); % Equivalent Compliance
    
    f0=1/(2*pi)*sqrt((N^2+C0/CE)/(M*C0)); % Resonance Frequency of Uncoupled PZT Ring
    ka0 = ((2*pi.*f0)./sp).*Geometry.a; % ka of Uncoupled PZT Ring
    
    f0Save(:,gpnum) = f0;
    ka0Save(:,gpnum) = ka0;
    
    %% Find the Lists of Zeros of J_0'(x)=0
    % Checked accuracy until n=10000
    
    J0dot=@(x) -besselj(1,x);
    rootj0dot=zeros(1,length(f));
    
    for i=1:n;
        rootj0dot(i)=fzero(J0dot,[i i+1]*pi); % Need to ignore first zero x=0
    end
    
    disp('Calculating Material Parameters is Finished');

%% State variables and T Matrix of the system
[StateVector, TMatrix]=StateVaribles_TMatrix_Subrutine(zlength1,zlength2,num);
disp('Calculating State variables and T Matrix of the system are Finished');
    
    %% Circuit Parameters of Inner Fluid
    [z_rr,z_rz,z_1,z_2,z_PRm]=CavityImpedance(LaRatio,Geometry.t,Geometry.a,ka,rhopzt,rhowater,sp,sE11,n,rootj0dot);
    
    %% Circuit Parameters about Gap Fluid 
    Cavityimp.gz_rr = zeros(length(ka),num-1);
    Cavityimp.gz_rz = zeros(length(ka),num-1);
    Cavityimp.gz_1 = zeros(length(ka),num-1);
    Cavityimp.gz_2 = zeros(length(ka),num-1);

    for Count = 1 : num - 1
        [gz_rr,gz_rz,gz_1,gz_2,~]=CavityImpedance(gaRatio(Count),0,0,ka,rhopzt,rhowater,sp,sE11,n,rootj0dot);
        
        Cavityimp.gz_rr(:,Count) = gz_rr.' .* (Geometry.g(Count)/L);
        Cavityimp.gz_rz(:,Count) = gz_rz.' .* (Geometry.g(Count)/L);
        Cavityimp.gz_1(:,Count) = gz_1.' .* (Geometry.g(Count)/L);
        Cavityimp.gz_2(:,Count) = gz_2.' .* (Geometry.g(Count)/L);
    end  

disp('Calculating Circuit Impedance Parameters is Finished');

    %% Circuit Parameters about Radiation Impedance (For 2Array)
    
%     [z_RMatrix, exceldata]=RadiationImp2Array(rhowater,sp,ka,a,L,zlength);
%     disp('Importing Raiation Impedance Parameters is Finished');
    
    %% Circuit Parameters about Radiation Impedance (Table)
        
%     [z_RMatrix]=RadiationImp_ParameterStudy(rhowater,sp,ka,a,L,g,zlength,RadiationTable,TableParameter);
[z_RAMatrix, z_RMatrix]=Radiation_Impedance(rhowater,sp,f,Geometry,num);
disp('Calculating Raiation Impedance Parameters is Finished');
    
    %% Calculate Total Admittance
    
    clear Y;
    
    % Construct Matrix and Vector
    z_cMatrix = cell(1,length(ka)); % Cavity mode impedance matrix
    z_PRMatrix = cell(1,length(ka)); % Piezoelectric ring matrix
    for i=1:length(ka);
        
        z_cMatrix{1,i}=zeros(zlength,zlength);
        for j=1:zlength;
            
            % input diagonal data
            if (j==1)||(j==zlength); % First, end diagonal -> Z1
                z_cMatrix{1,i}(j,j)=z_1(i);
            elseif (j==3)||(j==5)||(j==11)||(j==13); % Odd diagonal except First, end -> Z1+Z1gap
%                 z_cMatrix{1,i}(j,j)=z_1(i)+z_1gap(i); 
                  z_cMatrix{1,i}(j,j)=z_1(i)+Cavityimp.gz_rr(i,1);
            elseif (j==7)||(j==9); 
                  z_cMatrix{1,i}(j,j)=z_1(i)+Cavityimp.gz_rr(i,2); 
            elseif mod(j,4)==2; % Radial(cavity)
                z_cMatrix{1,i}(j,j)=z_rr(i);  
            elseif (j==4)||(j==12); % Radial(Side gap)
                z_cMatrix{1,i}(j,j)=Cavityimp.gz_rr(i,1);
            elseif j==8; % Radial(Center gap)
                z_cMatrix{1,i}(j,j)=Cavityimp.gz_rr(i,2);
            end
            
            % input tridiagonal data
            if mod(j,4)==2; % data of cavity (Center : Zrr)
                z_cMatrix{1,i}(j-1,j)=-z_rz(i);
                z_cMatrix{1,i}(j,j-1)=-z_rz(i);
                z_cMatrix{1,i}(j-1,j+1)=-z_2(i);
                z_cMatrix{1,i}(j+1,j-1)=-z_2(i);
                z_cMatrix{1,i}(j,j+1)=z_rz(i);
                z_cMatrix{1,i}(j+1,j)=z_rz(i);
            elseif (j==4)||(j==12); % data of gap (Side : Zrrgap)
                z_cMatrix{1,i}(j-1,j)=-Cavityimp.gz_rz(i,1);
                z_cMatrix{1,i}(j,j-1)=-Cavityimp.gz_rz(i,1);
                z_cMatrix{1,i}(j-1,j+1)=-Cavityimp.gz_2(i,1);
                z_cMatrix{1,i}(j+1,j-1)=-Cavityimp.gz_2(i,1);
                z_cMatrix{1,i}(j,j+1)=Cavityimp.gz_rz(i,1);
                z_cMatrix{1,i}(j+1,j)=Cavityimp.gz_rz(i,1);
            elseif (j==8); % data of gap (Center : Zrrgap)
                z_cMatrix{1,i}(j-1,j)=-Cavityimp.gz_rz(i,2);
                z_cMatrix{1,i}(j,j-1)=-Cavityimp.gz_rz(i,2);
                z_cMatrix{1,i}(j-1,j+1)=-Cavityimp.gz_2(i,2);
                z_cMatrix{1,i}(j+1,j-1)=-Cavityimp.gz_2(i,2);
                z_cMatrix{1,i}(j,j+1)=Cavityimp.gz_rz(i,2);
                z_cMatrix{1,i}(j+1,j)=Cavityimp.gz_rz(i,2);
            end
        end
        
        z_PRMatrix{1,i}=zeros(zlength,zlength);
        for j=1:num;
            z_PRMatrix{1,i}(4*j-2,4*j-2)=z_PRm(i);
        end
    end
    
    VectorT2=zeros(zlength,num);
    for i=1:num
        VectorT2(4*i-2,i)=1;
    end
    
    Y = cell(1,num);
    Admittance_Factor = N^2/(rhowater*sp*2*pi*Geometry.a*L);
    InversMatrix = cell(1,length(ka));
    
    for i=1:num
        Y{1,i} = zeros(length(ka),1); %초기화
    end
    
    for i=1:length(ka);
        InversMatrix{1,i} = (inv(z_cMatrix{1,i} + z_PRMatrix{1,i} + z_RAMatrix{1,i}));
        Y_temp = (1i*(2*pi*f(i))*C0) + Admittance_Factor * (VectorT2' * InversMatrix{1,i} * VectorT2 * ones(num, 1));
        for j=1:num
            Y{1,j}(i,1) =  Y_temp(j,1);
        end
    end
    
    % Admittance Data Cell to double Matrix
    DoubleY(:,1) = Y{1,1}(:,1);
    DoubleY(:,2) = Y{1,2}(:,1);
    clear Y;
    Y = DoubleY;
    clear DoubleY;
    % Admittance Data Cell to double Matrix
    
    % NAN 찾아서 채우기
    % Ring 1
    CheckY1 = Y(:,1);
    CheckNANVector = isnan(CheckY1);
    [NANn,NANm]=find(CheckNANVector == 1);
    
    if isempty(NANm) == 0
        CheckY1(NANn) = interp1(ka,CheckY1,ka(NANn),'spline');
    end
    
    % Ring 2
    CheckY2 = Y(:,2);
    CheckNANVector = isnan(CheckY2);
    [NANn,NANm]=find(CheckNANVector == 1);
    
    if isempty(NANm) == 0
        CheckY2(NANn) = interp1(ka,CheckY2,ka(NANn),'spline');
    end
    
    Y(:,1) = CheckY1;
    Y(:,2) = CheckY2;
    
    clear CheckY1 CheckY2;
    % NAN 찾아서 채우기
    
    % gValue == 0.20 & L = 0.1925 계산 오류 보정
    if (Geometry.g(1,1) == 0.20) && (L == 0.1925)
        ErrNum = find(f == 3750);
        disp('NAN find!')
        Y(ErrNum,1) = nan;
        Y(ErrNum,1) = interp1(f,Y(:,1),f(ErrNum),'spline');
        Y(ErrNum,2) = nan;
        Y(ErrNum,2) = interp1(f,Y(:,2),f(ErrNum),'spline');
    end
    % gValue == 0.20 & L = 0.1925 계산 오류 보정
    
    % gValue == 0.24 & L = 0.1925 계산 오류 보정
    if (Geometry.g(1,1) == 0.24) && (L == 0.1925)
        ErrNum = find(f == 3125);
        disp('NAN find!')
        Y(ErrNum,1) = nan;
        Y(ErrNum,1) = interp1(f,Y(:,1),f(ErrNum),'spline');
        Y(ErrNum,2) = nan;
        Y(ErrNum,2) = interp1(f,Y(:,2),f(ErrNum),'spline');
    end
    % gValue == 0.24 & L = 0.1925 계산 오류 보정
    
    Y1Save(:,gpnum) = Y(:,1);
    Y2Save(:,gpnum) = Y(:,2);
    RealY = real(Y1Save(:,gpnum));
    
    figure(1)
    plot(f,RealY(:,1))
    ExportCircuit = [f', RealY(:,1)];
    grid on;
    %     StrL = num2str(L);
    title('Conductance','fontsize',20, 'fontangle','italic');
    xlabel('Frequency (Hz)','fontsize',20, 'fontangle','italic');
    ylabel('Conductance','fontsize',20, 'fontangle','italic');
    set(gca, 'fontsize',16)
    set(gcf, 'color', 'w')
    hold on
    
    % ka normalized plot
    KaRatio = ka./ka0;
    KaRatioSave(:,gpnum) = KaRatio';
    
    figure(11)
    plot(KaRatio, RealY(:,1))
    grid on;
    title('Conductance','fontsize',20, 'fontangle','italic');
    xlabel({'ka normalized by (ka)_r ','= [ (ka) / (ka)_r ]     '},'fontsize',20, 'fontangle','italic');
    ylabel('Conductance','fontsize',20, 'fontangle','italic');
    set(gca, 'fontsize',16)
    set(gcf, 'color', 'w')
    hold on
    
    pause(0.1)
    
    disp('Admittance Caclulation is Finished');
    
    % Temp Velocity
    u = cell(1,num);
    for i=1:length(ka);
        u_temp = InversMatrix{1,i} * VectorT2;
        for j=1:num
            u{1,j}(i,1) =  u_temp(4*j-2,j);
        end
    end

    %% Calculate TVR at (r,z) and 1V input
    disp('TVR Calculate is Started');

% Consturct Matrix Area
Matrix.Area = zeros(zlength2, zlength2);
RLCount = 1;
gCount = 1;
for Count = 1 : zlength2
    if Count == 1
        Matrix.Area(Count, Count) = pi*Geometry.a^2; % 필요 시 부호 수정
    elseif Count == zlength2
        Matrix.Area(Count, Count) = -pi*Geometry.a^2; % 필요 시 부호 수정
    elseif rem(Count,2) == 0
        Matrix.Area(Count, Count) = 2*pi*Geometry.a*Geometry.RL(RLCount,1);
        RLCount = RLCount+1;
    else
        Matrix.Area(Count, Count) = 2*pi*Geometry.a*Geometry.g(gCount,1);
        gCount = gCount+1;
    end
end
clear RLCount gCount Count
% Consturct Matrix Area

% Consturct Matrix Aspect Ratio
% Matrix.AR = (2*pi*Geometry.a*Geometry.L) * inv(Matrix.Area);
Matrix.AR = (2*pi*Geometry.a*Geometry.RL(1,1)) * inv(Matrix.Area);
    
 
    % Construct Matrix and Vector
    p_uMatrix = cell(1,length(ka)); % p_u Matrix
    p_pMatrix = cell(1,length(ka)); % p_wp Matrix
    
    for i=1:length(ka);        
        [p_uMatrix{1,i},p_pMatrix{1,i}]=P_parameter_Subrutine(R0,Z0,Geometry.Zam,ka(1,i),Geometry.a,Geometry.L,zlength2,Matrix.AR,num);        
    end
    
    TVR_p = zeros(length(ka),2);
    % InversMatrix= cell(1,length(ka));
    TVR_Factor = N/(2*pi*Geometry.a*L);
    
    for i=1:length(ka);
       
%         TVR_ptemp = (p_uMatrix{1,i}.' * eye(length(p_uMatrix{1,i})) + p_wpMatrix{1,i}.'* Matrix.AR * z_RMatrix{1,i})...
%             *TMatrix.Tsa* InversMatrix{1,i} * VectorT2 * ones(num, 1);
        TVR_ptemp = (p_uMatrix{1,i}.' * eye(length(p_uMatrix{1,i})) + p_pMatrix{1,i}.' * z_RAMatrix{1,i})...
            * InversMatrix{1,i} * VectorT2 * ones(num, 1);

        TVR_p(i,:) = TVR_Factor .* TVR_ptemp ;
    end
    
    TVR = zeros(length(ka),1);
    for i=1:length(ka);
        TVR(i,1)=20*log10(sqrt(R0^2+Z0^2)*sqrt(TVR_p(i,1)*conj(TVR_p(i,1))/2)/1e-6);
    end
    
    % NAN 찾아서 채우기
    CheckNANVector = isnan(TVR);
    [NANn,NANm]=find(CheckNANVector == 1);
    
    if isempty(NANm) == 0
        TVR(NANn) = interp1(ka,TVR,ka(NANn),'spline');
    end
    % NAN 찾아서 채우기
    
    % gValue == 0.20 & L = 0.1925 계산 오류 보정
    if (Geometry.g(1,1) == 0.20) && (L == 0.1925)
        ErrNum = find(f == 3750);
        disp('NAN find!')
        TVR(ErrNum) = nan;
        TVR(ErrNum) = interp1(f,TVR,f(ErrNum),'spline');
    end
    % gValue == 0.20 & L = 0.1925 계산 오류 보정
    
    % gValue == 0.24 & L = 0.1925 계산 오류 보정
    if (Geometry.g(1,1) == 0.24) && (L == 0.1925)
        ErrNum = find(f == 3125);
        disp('NAN find!')
        TVR(ErrNum) = nan;
        TVR(ErrNum) = interp1(f,TVR,f(ErrNum),'spline');
    end
    % gValue == 0.24 & L = 0.1925 계산 오류 보정
    
    TVRSave(:,gpnum) = TVR.';
    
    figure(2)
    plot(f,TVR);
    grid on;
    title('TVR','fontsize',20, 'fontangle','italic');
    xlabel('Frequency (Hz)','fontsize',20, 'fontangle','italic');
    ylabel('TVR(dB)','fontsize',20, 'fontangle','italic');
    set(gca, 'fontsize',16)
    set(gcf, 'color', 'w')
    hold on
    
    figure(22)
    plot(KaRatio, TVR)
    grid on;
    title('TVR','fontsize',20, 'fontangle','italic');
    xlabel({'ka normalized by (ka)_r ','= [ (ka) / (ka)_r ]     '},'fontsize',20, 'fontangle','italic');
    ylabel('TVR(dB)','fontsize',20, 'fontangle','italic');
    set(gca, 'fontsize',16)
    set(gcf, 'color', 'w')
    hold on
        
    disp('TVR Caculation is Finished');
    pause(0.1)
    


end
hold off
toc
%% Contourf Plot
% Interpolation 으로 그림 resolution 증가하기
% SurPlotLaRatioSpline = L_ParameterSpline/a;

figure(3)
% SurPlotgaRatio = g_Parameter_Side/Geometry.a;
% SurPlotgLRatio = g_Parameter_Side/L;
SurPlotgaRatio = g_Parameter_Center/Geometry.a;
SurPlotgLRatio = g_Parameter_Center/L;
[MeshSurPLotka0, MeshSurPlotga] = meshgrid(KaRatio,SurPlotgaRatio);
[MeshSurPLotka1, MeshSurPlotgL] = meshgrid(KaRatio,SurPlotgLRatio);
SurRealYData = real(Y1Save.');
SurImagYData = imag(Y1Save.');

% SurPlotLaRatioSpline = 1000;
% [MeshSurPLotka0Spline, MeshSurPlotLaSpline] = meshgrid(KaRatio,SurPlotLaRatioSpline);
% SurRealYDataSpline = interp2(MeshSurPLotka0,MeshSurPlotLa,SurRealYData,MeshSurPLotka0Spline, MeshSurPlotLaSpline,'spline');

contourf(MeshSurPLotka0,MeshSurPlotga,SurRealYData,80 ,'LineStyle', 'none' ), colorbar
% contourf(MeshSurPLotka0Spline,MeshSurPlotLaSpline,SurRealYDataSpline,'LineStyle', 'none' ), colorbar
title('Conductance','fontsize',20);
xlabel({'(ka) / (k_0a) '},'fontsize',20, 'fontangle','italic');
ylabel('Gap Ratio (g/a)','fontsize',20, 'fontangle','italic');
set(gca, 'fontsize',16)
set(gcf, 'color', 'w')
grid off
% 
% figure(33)
% mesh(MeshSurPLotka0,MeshSurPlotga,SurRealYData), colorbar
% % mesh(MeshSurPLotka0Spline,MeshSurPlotLaSpline,SurRealYDataSpline), colorbar
% zlabel('Conductance','fontsize',20, 'fontangle','italic');
% xlabel({'ka normalized by (ka)_0 ','= [ (ka) / (ka)_0 ]     '},'fontsize',20, 'fontangle','italic');
% ylabel('Aspect Ratio (L/a)','fontsize',20, 'fontangle','italic');
% set(gca, 'fontsize',16)
% set(gcf, 'color', 'w')

figure(333)
surf(MeshSurPLotka0,MeshSurPlotga,SurRealYData), colorbar
% mesh(MeshSurPLotka0Spline,MeshSurPlotLaSpline,SurTVRDataSpline), colorbar
zlabel('Conductance','fontsize',20, 'fontangle','italic');
xlabel({'ka normalized by (ka)_0 ','= [ (ka) / (ka)_0 ]     '},'fontsize',20, 'fontangle','italic');
ylabel('Ratio (g/a)','fontsize',20, 'fontangle','italic');
set(gca, 'fontsize',16)
set(gcf, 'color', 'w')
shading interp

figure(334)
surf(MeshSurPLotka1,MeshSurPlotgL,SurRealYData), colorbar
% mesh(MeshSurPLotka0Spline,MeshSurPlotLaSpline,SurTVRDataSpline), colorbar
zlabel('Conductance','fontsize',20, 'fontangle','italic');
xlabel({'ka normalized by (ka)_0 ','= [ (ka) / (ka)_0 ]     '},'fontsize',20, 'fontangle','italic');
ylabel('Ratio (g/L)','fontsize',20, 'fontangle','italic');
set(gca, 'fontsize',16)
set(gcf, 'color', 'w')
shading interp

% figure(334)
% surf(MeshSurPLotka0,MeshSurPlotga,SurImagYData), colorbar
% % mesh(MeshSurPLotka0Spline,MeshSurPlotLaSpline,SurTVRDataSpline), colorbar
% zlabel('Susceptance','fontsize',20, 'fontangle','italic');
% xlabel({'ka normalized by (ka)_0 ','= [ (ka) / (ka)_0 ]     '},'fontsize',20, 'fontangle','italic');
% ylabel('Aspect Ratio (L/a)','fontsize',20, 'fontangle','italic');
% set(gca, 'fontsize',16)
% set(gcf, 'color', 'w')
% shading interp

% figure(4)
SurTVRData = TVRSave';

% SurTVRDataSpline = interp2(MeshSurPLotka0,MeshSurPlotLa,SurTVRData,MeshSurPLotka0Spline, MeshSurPlotLaSpline,'spline');

% contourf(MeshSurPLotka0,MeshSurPlotga,SurTVRData,80), colorbar
% caxis([80,140])
% % contourf(MeshSurPLotka0Spline,MeshSurPlotLaSpline,SurTVRDataSpline), colorbar
% title('TVR [dB]','fontsize',20, 'fontangle','italic');
% xlabel({'ka normalized by (ka)_0 ','= [ (ka) / (ka)_0 ]     '},'fontsize',20, 'fontangle','italic');
% ylabel('Aspect Ratio (L/a)','fontsize',20, 'fontangle','italic');
% set(gca, 'fontsize',16)
% set(gcf, 'color', 'w')
% grid on

% figure(44)
% mesh(MeshSurPLotka0,MeshSurPlotLa,SurTVRData), colorbar
% % mesh(MeshSurPLotka0Spline,MeshSurPlotLaSpline,SurTVRDataSpline), colorbar
% zlabel('TVR [dB]','fontsize',20, 'fontangle','italic');
% xlabel({'ka normalized by (ka)_0 ','= [ (ka) / (ka)_0 ]     '},'fontsize',20, 'fontangle','italic');
% ylabel('Aspect Ratio (L/a)','fontsize',20, 'fontangle','italic');
% set(gca, 'fontsize',16)
% set(gcf, 'color', 'w')

figure(444)
surf(MeshSurPLotka0,MeshSurPlotga,SurTVRData), colorbar
caxis([80,140])
% mesh(MeshSurPLotka0Spline,MeshSurPlotLaSpline,SurTVRDataSpline), colorbar
zlabel('TVR [dB]','fontsize',20);
xlabel('(ka) / (k_0a)','fontsize',20, 'fontangle','italic');
ylabel('Gap Ratio (g/a)','fontsize',20, 'fontangle','italic');
set(gca, 'fontsize',16)
set(gcf, 'color', 'w')
shading interp

figure(445)
surf(MeshSurPLotka1,MeshSurPlotgL,SurTVRData), colorbar
caxis([80,140])
% mesh(MeshSurPLotka0Spline,MeshSurPlotLaSpline,SurTVRDataSpline), colorbar
zlabel('TVR [dB]','fontsize',20, 'fontangle','italic');
xlabel({'ka normalized by (ka)_0 ','= [ (ka) / (ka)_0 ]     '},'fontsize',20, 'fontangle','italic');
ylabel('Ratio (g/L)','fontsize',20, 'fontangle','italic');
set(gca, 'fontsize',16)
set(gcf, 'color', 'w')
shading interp

%% Calculate Directivity of R=RD & f=fD at x-Axis
% 
% % Define sets of the points to calculate
% RD=input('Radius to Measure (m) = '); % Radius to measure
% fD=input('Frequency to Measure (Hz) = '); % Frequency to measure
% i=1;
% while(i<=length(f))  % Measured frequency should be the element of former f data
%    if f(i)==fD;
%        break;
%    end
%    i=i+1;
%    if i>length(f)
%        disp('Input the valid frequency value');
%        fD=input('Frequency to Measure (Hz) = '); % Frequency to measure
%        i=1;
%    end
% end
% % numD=input('Number of Measuring Points (#) = '); % Number of points
% numD=181;
% 
% % Define sets of the points to calculate
% RDr=zeros(1,length(numD)); % r-Points Data
% RDz=zeros(1,length(numD)); % z-Points Data
% 
% theta=linspace(-pi/2,pi/2,numD); % -pi/2 to pi/2
% thetadeg=theta*180/pi; % Angle (degree)
% for i=1:numD; % Points
%     RDr(i)=RD*cos(theta(i));
%     RDz(i)=RD*sin(theta(i));
% end
% 
% % find the points fi when f(fi)=fD
% for i=1:length(f);
%     if f(i)==fD;
%         fi=i;
%     end
% end
% 
% Dka = ((2*pi*fD)/sp)*a; % ka, dimensionless parameter by Been
% TVR_Factor = N/(2*pi*a*L);
% Pdata=zeros(1,length(theta)); % Total Pressure Data
% TVRdata=zeros(1,length(theta)); % TVR Datas
% 
% for iD=1:numD
%     [Dp_1u, Dp_r1u, Dp_r2u, Dp_2u, Dp_gu, Dp_1p, Dp_r1p, Dp_r2p, Dp_gp, Dp_2p]=Pres_para_2Array(RDr(iD),RDz(iD),Dka,LaRatio,a,L,g);
%     % Construct Matrix and Vector
%         
%     Dp_uMatrix = [Dp_1u; Dp_r1u; 0; Dp_gu; 0; Dp_r2u; Dp_2u];
%     Dp_wpMatrix = [-(2*LaRatio)*Dp_1p; Dp_r1p; 0; (g/L).* Dp_gp; 0; Dp_r2p; (2*LaRatio)*Dp_2p]; % V2
%     
%     DTVR_ptemp = (Dp_uMatrix.' * eye(length(Dp_uMatrix)) + Dp_wpMatrix.' * z_RMatrix{1,fi})...
%         * InversMatrix{1,fi} * VectorT2 * [1;1];
%     Pdata(iD) = TVR_Factor .* DTVR_ptemp ;
%         
%     TVRdata(iD)=20*log10(sqrt(RDr(iD)^2+RDz(iD)^2)*sqrt(Pdata(iD)*conj(Pdata(iD))/2)/1e-6);
% end
% 
% % Pressure of axis (for calculate beam pattern)
% [Dp_1u, Dp_r1u, Dp_r2u, Dp_2u, Dp_gu, Dp_1p, Dp_r1p, Dp_r2p, Dp_gp, Dp_2p]=Pres_para_2Array(RD,0,Dka,LaRatio,a,L,g);
% 
% Dp_uMatrix = [Dp_1u; Dp_r1u; 0; Dp_gu; 0; Dp_r2u; Dp_2u];
% Dp_wpMatrix = [-(2*LaRatio)*Dp_1p; Dp_r1p; 0; (g/L).* Dp_gp; 0; Dp_r2p; (2*LaRatio)*Dp_2p]; % V2
% 
% DTVR_ptemp = (Dp_uMatrix.' * eye(length(Dp_uMatrix)) + Dp_wpMatrix.' * z_RMatrix{1,fi})...
%     * InversMatrix{1,fi} * VectorT2 * [1;1];
% Paxisdata = TVR_Factor .* DTVR_ptemp ;
% 
% Ddata=20*log10(abs(Pdata)./abs(Paxisdata)); % Calculate Beam Pattern
% figure (3)
% dirplot(thetadeg,Ddata); % Plot polar plot based on the plotting programs
% ExpertD = [thetadeg',Ddata'];

%% Radiation Resistance

% for i = 1:31
%     j = 20*i+41;
%     Z_11(i,1) = z_RAMatrix{1,j}(1,1);
%     Z_12(i,1) = z_RAMatrix{1,j}(1,2);
%     Z_13(i,1) = z_RAMatrix{1,j}(1,4);
%     Z_14(i,1) = z_RAMatrix{1,j}(1,6);
%     Z_15(i,1) = z_RAMatrix{1,j}(1,7);
%     Z_22(i,1) = z_RAMatrix{1,j}(2,2);
%     Z_23(i,1) = z_RAMatrix{1,j}(2,4);
%     Z_24(i,1) = z_RAMatrix{1,j}(2,6);
%     Z_33(i,1) = z_RAMatrix{1,j}(4,4);
% end
% Z_temp = [Z_11 Z_12 Z_13 Z_14 Z_15 Z_22 Z_23 Z_24 Z_33];
% Z = [real(Z_temp) imag(Z_temp)];
