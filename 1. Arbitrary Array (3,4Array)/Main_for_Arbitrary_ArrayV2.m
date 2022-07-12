%% FFR Circuit Model Code (Arbitrary Array)
% Made by Kyounghun Been

clc
clear
close all

tic
%% Input Geometry and Material Inputs of FFR
% Assume Matrial of PZT is PZT-5H and Medium is Water
% PZT Transducetion Parameters (PZT-5H)
Material.sE11=1.65e-11;
Material.eT33=3400*8.8541878176*10^-12;
Material.d31=-2.74e-10;
Material.k31=sqrt(Material.d31^2/(Material.sE11*Material.eT33));
Material.eS33=Material.eT33*(1-Material.k31^2);

Material.rhopzt=7500; % Density of PZT
Material.rhowater=999; % Density of Water
Material.sp=1500; % Sound Speed at Water

NumR=input('Number of FFR = '); % Number of FFR
if NumR == 1
    Wave.f=25:1:4000; % Set Frequency Range single
else
    Wave.f=25:5:4000;  % Set Frequency Range for Calculate 2Array
end

% f=25:5:4000;  % Set Frequency Range for Calculate 2Array
% f=25:1:4000; % Set Frequency Range single
Wave.omega=2*pi*Wave.f; % Angular Frequency
% k=@(w) w/sp; % Wave Number

% Input the coordinate to mesure
% R0=input('r-Distance (m) = '); % the r-distance to calculate TVR
% Z0=input('z-Distance (m) = '); % the z-distance to calculate TVR
% while (R0<=a)&&(abs(Z0)<=(N*L+(N-1)*g)/2); % The point should be at the outside of the Array
%     disp('Input valid value of R which is R>a');
%     R0=input('r-Distance (m) = ');
%     Z0=input('z-Distance (m) = ');
% end

% Input the coordinate to mesure
TVR.R0=2; % the r-distance to calculate TVR
TVR.Z0=0; % the z-distance to calculate TVR

% Input the Voltage
Voltage = zeros(NumR,1); % Voltage
for Count = 1 : NumR
    Voltage(Count,1) = 1;
end

zlength1=4*NumR-1; % length of the matrix of array number = N
zlength2=2*NumR+1; % length of the matrix of array number = N

Geometry.a=0.175; % Ring Average Radius (Fix)
Geometry.t=0.04; % Ring Thickness

Geometry.RL = zeros(NumR,1); % Ring Height
for Count = 1 : NumR
    Geometry.RL(Count,1) = 0.1925;
end
clear Count

if NumR==1; % Single FFR
    Geometry.g=0;
else % Array FFR
    Geometry.g = zeros(NumR-1,1); % Gab Height
    for Count = 1 : NumR - 1
        Geometry.g(Count,1) = 0.08;
    end
    clear Count
end

Geometry.L = sum(Geometry.RL) + sum(Geometry.g); % Total length of FFR array 
Wave.ka = ((2*pi.*Wave.f)./Material.sp).*Geometry.a; % ka, dimensionless parameter by Been
Geometry.LaRatio = Geometry.RL./Geometry.a; % L/a ratio, dimensionless parameter by Been
Geometry.gaRatio = Geometry.g./Geometry.a; % g/a ratio, dimensionless parameter by Been
% Dimfactor = Material.rhowater*Material.sp*2*pi*a.*L;
% DimfactorGap = Material.rhowater*Material.sp*2*pi*a.*g;

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

%% Calculation Admittance Factor and p0 Factor
Factor.Admittance = Material.rhowater*Material.sp*2*pi*Geometry.a*Geometry.L;
Factor.p0 = 1/(2*pi*Geometry.a*Geometry.L);
% Factor.p0 = Material.d31 / (Geometry.a * Material.sE11);

%% Circuit Parameters of PZT
% Ignore Electrical Dissipation G0 and Mechanical Damping Rm

PZT.C0=(2*pi*Geometry.a.*Geometry.RL./Geometry.t).*Material.eS33; % Equivalent Clamped Capacitance
PZT.N=(2*pi.*Geometry.RL).*(Material.d31/Material.sE11); % Equivalent Electromechanical turns Ratio
PZT.M=Material.rhopzt*2*pi*Geometry.a*Geometry.RL(1)*Geometry.t; % Mass of PZT
PZT.CE=(Material.sE11*Geometry.a)/(2*pi*Geometry.t*Geometry.RL(1)); % Equivalent Compliance

Wave.f0=1/(2*pi).*sqrt((PZT.N.^2+PZT.C0./PZT.CE)./(PZT.M.*PZT.C0)); % Resonance Frequency of Uncoupled PZT Ring
Wave.ka0 = ((2*pi.*Wave.f0)./Material.sp).*Geometry.a; % ka of Uncoupled PZT Ring

Wave.KaRatio = zeros(length(Wave.ka),NumR);
for Count=1:NumR
    Wave.KaRatio(:,Count) = Wave.ka'./Wave.ka0(Count,1); % For nomalized ka0
end
clear Count

%% Find the Lists of Zeros of J_0'(x)=0
% Checked accuracy until n=10000

J0dot=@(x) -besselj(1,x);
rootj0dot=zeros(1,length(Wave.f));

for Count=1:n;
    rootj0dot(Count)=fzero(J0dot,[Count Count+1]*pi); % Need to ignore first zero x=0
end
clear Count

disp('Calculating Material Parameters is Finished');

%% State variables and T Matrix of the system
[StateVector, TMatrix]=StateVaribles_TMatrix_Subrutine(zlength1,zlength2,NumR);
disp('Calculating State variables and T Matrix of the system are Finished');

%% Circuit Parameters about Radiation Impedance
if NumR == 1
    [z_RMatrix,~]=Radiation_Impedance_Single_Subrutine(Material.rhowater,Material.sp,Wave.ka,Geometry.a,Geometry.L,zlength1);
elseif NumR == 2
    [z_RMatrix, ~]=Radiation_Impedance_2Array_Subrutine(Material.rhowater,Material.sp,Wave.ka,Geometry.a,Geometry.L,zlength1);
elseif NumR == 3
    [z_RMatrix, ~]=Radiation_Impedance_3Array_Subrutine(Material.rhowater,Material.sp,Wave.ka,Geometry.a,Geometry.L,zlength1);
elseif NumR == 4
    [z_RMatrix, ~]=Radiation_Impedance_4Array_Subrutine(Material.rhowater,Material.sp,Wave.ka,Geometry.a,Geometry.L,zlength1);
end
disp('Importing Raiation Impedance Parameters is Finished');

% %% Circuit Parameters about Radiation Impedance (For Single)
% [z_RMatrix,exceldata]=RadiationImpedance(rhowater,sp,ka,a,L,zlength1);
% disp('Importing Raiation Impedance Parameters is Finished');

% %% Circuit Parameters about Radiation Impedance (For 2Array)
% [z_RMatrix, exceldata]=RadiationImp2Array(rhowater,sp,ka,a,L,zlength1);
% disp('Importing Raiation Impedance Parameters is Finished');

%% Circuit Parameters of Inner Fluid
Matrix.zRc = zeros(3,3,NumR);
Matrix.zgc = zeros(3,3,NumR-1);
Matrix.zPRm = zeros(3,3,NumR);

Cavityimp.Rz_rr = zeros(length(Wave.ka),NumR);
Cavityimp.Rz_rz = zeros(length(Wave.ka),NumR);
Cavityimp.Rz_1 = zeros(length(Wave.ka),NumR);
Cavityimp.Rz_2 = zeros(length(Wave.ka),NumR);
Cavityimp.Rz_PRm = zeros(length(Wave.ka),NumR);
Cavityimp.gz_rr = zeros(length(Wave.ka),NumR-1);
Cavityimp.gz_rz = zeros(length(Wave.ka),NumR-1);
Cavityimp.gz_1 = zeros(length(Wave.ka),NumR-1);
Cavityimp.gz_2 = zeros(length(Wave.ka),NumR-1);

for Count = 1 : NumR
    [z_rr,z_rz,z_1,z_2,z_PRm]=Cavity_Impedance_Subrutine(Geometry.LaRatio(Count),Geometry.t,Geometry.a,Wave.ka,Material.rhopzt,...
        Material.rhowater,Material.sp,Material.sE11,n,rootj0dot);
    
    Cavityimp.Rz_rr(:,Count) = z_rr.' .* (Geometry.RL(Count)/Geometry.L);
    Cavityimp.Rz_rz(:,Count) = z_rz.' .* (Geometry.RL(Count)/Geometry.L);
    Cavityimp.Rz_1(:,Count) = z_1.' .* (Geometry.RL(Count)/Geometry.L);
    Cavityimp.Rz_2(:,Count) = z_2.' .* (Geometry.RL(Count)/Geometry.L);
    Cavityimp.Rz_PRm(:,Count) = z_PRm.' .* (Geometry.RL(Count)/Geometry.L);
    
    % 비교를 위한 non-dim 변환
    Dimfactor = Material.rhowater*Material.sp*2*pi*Geometry.a*Geometry.L;
    Z_rr = z_rr.*Dimfactor;
    Z_rz = z_rz.*Dimfactor;
    Z_1 = z_1.*Dimfactor;
    Z_2 = z_2.*Dimfactor;
    Z_PRm = z_PRm.*Dimfactor;

end
% clear Count z_rr z_rz z_1 z_2 z_PRm
clear Count z_PRm

% Temp
ExportCavityImp = [imag(Cavityimp.Rz_rr), imag(Cavityimp.Rz_rz), imag(Cavityimp.Rz_1), imag(Cavityimp.Rz_2)]; 
% Temp

for Count = 1 : NumR - 1
    [z_rr,z_rz,z_1,z_2,~]=Cavity_Impedance_Subrutine(Geometry.gaRatio(Count),0,0,Wave.ka,Material.rhopzt,...
        Material.rhowater,Material.sp,Material.sE11,n,rootj0dot);
    
    Cavityimp.gz_rr(:,Count) = z_rr.' .* (Geometry.g(Count)/Geometry.L);
    Cavityimp.gz_rz(:,Count) = z_rz.' .* (Geometry.g(Count)/Geometry.L);
    Cavityimp.gz_1(:,Count) = z_1.' .* (Geometry.g(Count)/Geometry.L);
    Cavityimp.gz_2(:,Count) = z_2.' .* (Geometry.g(Count)/Geometry.L);
end
% clear Count z_rr z_rz z_1 z_2
clear Count

%% Constant Matrix Constuction
% Consturct Matrix TN
TNTerm2 = zeros(zlength1,zlength1);
for Count = 1 : NumR
    TNTerm1 = TMatrix.Tsre(:,:,Count)'*TMatrix.Tsre(:,:,Count);
    TNTerm2 = TNTerm2 + TNTerm1;
end
Matrix.TN = TNTerm2*TMatrix.Tc;
clear TNTerm1 TNTerm2 Count
% Consturct Matrix TN
    
% Consturct Matrix CO
Matrix.C0 = zeros(NumR, NumR);
for Count = 1 : NumR
    Matrix.C0(Count,Count) = PZT.C0(Count);
end
clear Count
% Consturct Matrix CO

% Consturct Matrix Nr
Matrix.Nr = zeros(NumR, NumR);
for Count = 1 : NumR
    Matrix.Nr(Count,Count) = PZT.N(Count);
end
clear Count
% Consturct Matrix Nr

% Consturct Matrix Voltage
Matrix.Voltage = ones(NumR, 1);
for Count = 1 : NumR
    Matrix.Voltage(Count,1) = Matrix.Voltage(Count,1) * Voltage(Count);
end
clear Count Voltage
% Consturct Matrix Voltage

%% Admittance
Admittance.Y = cell(1,NumR);
for Count = 1 : NumR
    Admittance.Y{1,Count} = zeros(length(Wave.ka),1); %초기화
end

disp('Admittance Caclulation is Started');
for MainCount = 1 : length(Wave.ka);
    
    % Construct Matrix zaR
    Matrix.zaR = z_RMatrix{1,MainCount};
    Matrix.zaRSave(:,:,MainCount) = Matrix.zaR;
    % Construct Matrix zaR 
    
    % Construct Matrix zR
    Matrix.zR = Matrix.zaR;
    MinusMatrix = 0; % 행이 감소함에 따라 삭제 해야 하는 행의 값을 보정해주는 변수
    for Count = 3 : zlength1-2
        if rem(Count,2) == 1
            Matrix.zR(Count-MinusMatrix,:) = [];
            Matrix.zR(:,Count-MinusMatrix) = [];
            MinusMatrix = MinusMatrix+1;
        end
    end
    Matrix.zRSave(:,:,MainCount) = Matrix.zR;
    clear Count MinusMatrix
    % Construct Matrix zR
    
    % Consturct Matrix zc
    for Count = 1 : NumR
        Matrix.zRc(:,:,Count) = ...
            [Cavityimp.Rz_1(MainCount,Count), Cavityimp.Rz_rz(MainCount,Count), -Cavityimp.Rz_2(MainCount,Count);...
            Cavityimp.Rz_rz(MainCount,Count), Cavityimp.Rz_rr(MainCount,Count), -Cavityimp.Rz_rz(MainCount,Count);...
            -Cavityimp.Rz_2(MainCount,Count), -Cavityimp.Rz_rz(MainCount,Count), Cavityimp.Rz_1(MainCount,Count)];
    end
    clear Count
    % Consturct Matrix zc
    
    % Consturct Matrix zPR
    for Count = 1 : NumR
        Matrix.zPRm(:,:,Count) = [0 0 0; 0 Cavityimp.Rz_PRm(MainCount,Count) 0; 0 0 0];
    end
    clear Count
    % Consturct Matrix zPR
    
    % Consturct Matrix zgc
    for Count = 1 : NumR-1
        Matrix.zgc(:,:,Count) = ...
            [Cavityimp.gz_1(MainCount,Count), Cavityimp.gz_rz(MainCount,Count), -Cavityimp.gz_2(MainCount,Count);...
            Cavityimp.gz_rz(MainCount,Count), Cavityimp.gz_rr(MainCount,Count), -Cavityimp.gz_rz(MainCount,Count);...
            -Cavityimp.gz_2(MainCount,Count), -Cavityimp.gz_rz(MainCount,Count), Cavityimp.gz_1(MainCount,Count)];
    end
    clear Count
    % Consturct Matrix zgc
    
    % Consturct Matrix za
    zaTerm1_temp = zeros(zlength1,zlength1,NumR);
    zaTerm1 = zeros(zlength1,zlength1);
    zaTerm2_temp = zeros(zlength1,zlength1,NumR-1);
    zaTerm2 = zeros(zlength1,zlength1);
    for Count = 1 : NumR
        zaTerm1_temp(:,:,Count) = TMatrix.Tsre(:,:,Count)'*...
            (Matrix.zRc(:,:,Count) + Matrix.zPRm(:,:,Count))*TMatrix.Tsre(:,:,Count);
        zaTerm1 = zaTerm1 + zaTerm1_temp(:,:,Count);
    end
    for Count = 1 : NumR - 1
        zaTerm2_temp(:,:,Count) = TMatrix.Tsge(:,:,Count)'*...
            Matrix.zgc(:,:,Count)*TMatrix.Tsge(:,:,Count);
        zaTerm2 = zaTerm2 + zaTerm2_temp(:,:,Count);
    end
    Matrix.za = zaTerm1 + zaTerm2;
    clear zaTerm1 zaTerm1_temp zaTerm2 zaTerm2_temp Count
    % Consturct Matrix za
    
    % Consturct Inverse Matrix [za + zr]
    Matrix.Inverse = inv(Matrix.za + Matrix.zaR);
    Matrix.InverseSave(:,:,MainCount) = Matrix.Inverse;
    % Consturct Inverse Matrix [za + zr]
            
    % Calculate Admittance    
    Admittance.Y_temp = ((Matrix.Nr' * Matrix.TN' * Matrix.Inverse * Matrix.TN * Matrix.Nr) / Factor.Admittance +...
        (1i*(2*pi*Wave.f(MainCount))*Matrix.C0)) * Matrix.Voltage;
    
    for Count=1:NumR
        Admittance.Y{1,Count}(MainCount,1) =  Admittance.Y_temp(Count,1);
    end
    clear Count
    % Calculate Admittance
end
% clear z_RMatrix
disp('Admittance Caclulation is Finished');

%% Admittance Plot
figure(1)
for Count=1:NumR
    plot(Wave.f,real(Admittance.Y{1,Count}))
    hold on
end
clear Count
hold off
grid on
title('Conductance','fontsize',20, 'fontangle','italic');
xlabel('Frequency [Hz]','fontsize',20, 'fontangle','italic');
ylabel('Conductance','fontsize',20, 'fontangle','italic');
set(gca, 'fontsize',16)
set(gcf, 'color', 'w')

% ka normalized plot
figure(11)
for Count=1:NumR
    plot(Wave.KaRatio(:,Count),real(Admittance.Y{1,Count}))
    hold on
end
clear Count
hold off
grid on
title('Conductance','fontsize',20, 'fontangle','italic');
xlabel({'ka normalized by (ka)_0 ','= [ (ka) / (ka)_0 ]     '},'fontsize',20, 'fontangle','italic');
ylabel('Conductance','fontsize',20, 'fontangle','italic');
set(gca, 'fontsize',16)
set(gcf, 'color', 'w')

%% Output
AdmittanceResult = [Wave.f', Wave.KaRatio, real(Admittance.Y{1,1})];

pause(0.1)

%% TVR
disp('TVR Caculation is Started');
% Consturct Matrix Area
Matrix.Area = zeros(zlength2, zlength2);
RLCount = 1;
gCount = 1;
for Count = 1 : zlength2
    if Count == 1
        Matrix.Area(Count, Count) = pi*Geometry.a^2;
    elseif Count == zlength2
        Matrix.Area(Count, Count) = -pi*Geometry.a^2;
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
Matrix.AR = (2*pi*Geometry.a*Geometry.L) * inv(Matrix.Area);
% Consturct Matrix Aspect Ratio

TVR.p0 = zeros(length(Wave.ka),1);
TVR.TVR = zeros(length(Wave.ka),1);
hwait = waitbar(0,'wait','Name','Calculating TVR');
for MainCount = 1 : length(Wave.ka);
    waitbar(MainCount/length(Wave.ka),hwait,sprintf('%1.0f %%',(MainCount/length(Wave.ka))*100))
    % Construct Matrix pu and pp
    Matrix.pu = zeros(zlength2,1);
    Matrix.pp = zeros(zlength2,1);
    for Count = 1 : zlength2
        [TVR.pu,TVR.pp]=P_parameter_Subrutine(TVR.R0,TVR.Z0,Geometry.Zam,Wave.ka(1,MainCount),Geometry.a,Geometry.L,Count,zlength2);
        Matrix.pu(Count,1) = TVR.pu;
        Matrix.pp(Count,1) = TVR.pp;
    end
    Matrix.puSave(:,:,MainCount) = Matrix.pu;
    Matrix.ppSave(:,:,MainCount) = Matrix.pp;
    % Construct Matrix pu and pp    
    
    TVR.p0_temp = (Matrix.pu.' * eye(length(Matrix.pu)) + Matrix.pp.' * Matrix.AR * Matrix.zRSave(:,:,MainCount))...
        * TMatrix.Tsa * Matrix.InverseSave(:,:,MainCount) * Matrix.TN * Matrix.Nr * Matrix.Voltage;
    TVR.p0(MainCount,:) = Factor.p0 .* TVR.p0_temp ;
        
    TVR.TVR(MainCount,1)=20*log10(sqrt(TVR.R0^2+TVR.Z0^2)*sqrt(TVR.p0(MainCount,1)*conj(TVR.p0(MainCount,1))/2)/1e-6);
end
delete(hwait)
clear MainCount Count
disp('TVR Caculation is Finished');

%% TVR Plot
figure(2)
plot(Wave.f,TVR.TVR(:,1));
grid on
title('TVR','fontsize',20, 'fontangle','italic');
xlabel('Frequency [Hz]','fontsize',20, 'fontangle','italic');
ylabel('TVR [dB]','fontsize',20, 'fontangle','italic');
set(gca, 'fontsize',16)
set(gcf, 'color', 'w')

% ka normalized plot
figure(22)
plot(Wave.KaRatio(:,1),TVR.TVR(:,1));
grid on
title('TVR','fontsize',20, 'fontangle','italic');
xlabel({'ka normalized by (ka)_0 ','= [ (ka) / (ka)_0 ]     '},'fontsize',20, 'fontangle','italic');
ylabel('TVR [dB]','fontsize',20, 'fontangle','italic');
set(gca, 'fontsize',16)
set(gcf, 'color', 'w')

toc