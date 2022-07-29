function [GD_ob, GW_ob] = HKI_Sub_ObPres_Mat(ro,zo,freq,c0,rho0,radius_a,length_l,NN)
%% Geometry & Boundary Condition Parameter
LaRatio = length_l/radius_a;

%% Define, number of node
Nb = NN(1); Nr = NN(2); Nt = NN(3);

TotNumPre = (Nb + Nr + Nt) - 2; % Total number of Pressure Vector
TotNumVel = (Nb + Nr + Nt); % Total number of Velocity Vector

detN_b = 1/(Nb-1); detN_r = 1/(Nr-1); detN_t = 1/(Nt-1); % gap of nodes

% 적분 구간 설정시 특이점 제외를 위한 Error 정도 설정, 0.1%
IntError_b = detN_b*0.002; IntError_r = detN_r*0.001; IntError_t = detN_t*0.002; 

%% Construct Node
% Node_eta_b = linspace(0,1,Nb)'; Node_eta_t = linspace(1,0,Nt)'; Node_zta_r = linspace(-0.5,0.5,Nr)';
Node_eta_b = (0:detN_b:1)'; Node_eta_t = (1:-detN_t:0)'; Node_zta_r = (-0.5:detN_r:0.5)'; 

%% Medium Parameter
% Water
c_water = c0; rho0_water = rho0;
ka = (2*pi*freq/c_water)*radius_a;

%% Position of observation
% ro = 1;
% zo = 0;

etaObNode = ro./radius_a;
ztaObNode = zo./length_l;
NumberOfObNode = length(etaObNode);

% Define Solid Angle
SolidAngle = 4*pi;

%% Cal Global D&W Matrix
% GD(PressureVectorNumber, VelocityVectorNumber)
GD_ob = zeros(NumberOfObNode,TotNumVel);
GW_ob = zeros(NumberOfObNode,TotNumPre);

MD_b = zeros(NumberOfObNode,Nb);
MD_r = zeros(NumberOfObNode,Nr);
MD_t = zeros(NumberOfObNode,Nt);

MW_b = zeros(NumberOfObNode,Nb);
MW_r = zeros(NumberOfObNode,Nr);
MW_t = zeros(NumberOfObNode,Nt);

parfor num = 1:NumberOfObNode
    % Define eto_o, zta_o
    eta_o = etaObNode(num);
    zta_o = ztaObNode(num);
    
    %% Integral Bottom, Surface 1, Global D
    R=@(eta,theta) HKI_Sub_R(eta_o,zta_o,eta,theta,-0.5,LaRatio); % Define R using Sub_R(eta_o,zta_o,eta,theta,zta,LaRatio)
    Int_D=@(eta,theta) HKI_Sub_IntFn_D(1,R,ka,eta,theta,0); % Define Integral Equation using Sub_IntFn_D(Surface Num,R,ka,eta,theta,zta)
    Fac=(-(1i*ka)/SolidAngle); % After Integral, multiply Fac.
    
    MD_b1 = zeros(1,Nb); % Temp Variable; MultiThread
    MD_b2 = zeros(1,Nb); % Temp Variable; MultiThread
        
    for num_b = 1:Nb
        IntRange = [Node_eta_b(num_b)-detN_b Node_eta_b(num_b)];
        Int_Fn=@(eta,theta) HKI_Sub_SFn(1,eta,IntRange).*Int_D(eta,theta); % Define Integral Fn, Sub_SFn is define Shape Function
        MD_b1(num_b) = HKI_Sub_Integral(Int_Fn,IntRange,IntError_b)*Fac;

        IntRange =  [Node_eta_b(num_b) Node_eta_b(num_b)+detN_b];
        Int_Fn=@(eta,theta) HKI_Sub_SFn(2,eta,IntRange).*Int_D(eta,theta); % Define Integral Fn
        MD_b2(num_b) = HKI_Sub_Integral(Int_Fn,IntRange,IntError_b)*Fac;
    end
    MD_b(num,:) = MD_b1 + MD_b2; % MD_b matrix construction; MultiThread
    
    %% Integral Ring, Surface 2, Global D
    R=@(zta,theta) HKI_Sub_R(eta_o,zta_o,1,theta,zta,LaRatio); % Define R
    Int_D=@(zta,theta) HKI_Sub_IntFn_D(2,R,ka,0,theta,zta); % Define Integral Equation
    Fac=((1i*ka)/SolidAngle)*LaRatio;
    
    MD_r1 = zeros(1,Nr); % Temp Variable; MultiThread
    MD_r2 = zeros(1,Nr); % Temp Variable; MultiThread
    
    for num_r = 1:Nr
        IntRange = [Node_zta_r(num_r)-detN_r Node_zta_r(num_r)];
        Int_Fn=@(zta,theta) HKI_Sub_SFn(1,zta+0.5,IntRange+0.5).*Int_D(zta,theta); % Define Integral Fn
        MD_r1(num_r) = HKI_Sub_Integral(Int_Fn,IntRange,IntError_r)*Fac;
        
        IntRange =  [Node_zta_r(num_r) Node_zta_r(num_r)+detN_r];
        Int_Fn=@(zta,theta) HKI_Sub_SFn(2,zta+0.5,IntRange+0.5).*Int_D(zta,theta); % Define Integral Fn
        MD_r2(num_r) = HKI_Sub_Integral(Int_Fn,IntRange,IntError_r)*Fac;
    end
    MD_r(num,:) = MD_r1 + MD_r2; % MD_b matrix construction; MultiThread
    
    %% Integral Top, Surface 3, Global D
    R=@(eta,theta) HKI_Sub_R(eta_o,zta_o,eta,theta,0.5,LaRatio); % Define R
    Int_D=@(eta,theta) HKI_Sub_IntFn_D(3,R,ka,eta,theta,0); % Define Integral Equation
    Fac=((1i*ka)/SolidAngle);
    
    MD_t1 = zeros(1,Nt); % Temp Variable; MultiThread
    MD_t2 = zeros(1,Nt); % Temp Variable; MultiThread
    
    for num_t = 1:Nt
        IntRange = [Node_eta_t(num_t)-detN_t Node_eta_t(num_t)];
        Int_Fn=@(eta,theta) HKI_Sub_SFn(1,eta,IntRange).*Int_D(eta,theta); % Define Integral Fn
        MD_t1(num_t) = HKI_Sub_Integral(Int_Fn,IntRange,IntError_t)*Fac;
        
        IntRange = [Node_eta_t(num_t) Node_eta_t(num_t)+detN_t];
        Int_Fn=@(eta,theta) HKI_Sub_SFn(2,eta,IntRange).*Int_D(eta,theta); % Define Integral Fn
        MD_t2(num_t) = HKI_Sub_Integral(Int_Fn,IntRange,IntError_t)*Fac;
    end
    MD_t(num,:) = MD_t1 + MD_t2; % MD_b matrix construction; MultiThread

    %% Integral Bottom, Surface 1, Global W
    R=@(eta,theta) HKI_Sub_R(eta_o,zta_o,eta,theta,-0.5,LaRatio); % Define R
    Int_W=@(eta,theta) HKI_Sub_IntFn_W(1,R,ka,0,eta,theta,0); % Define Integral Equation
    Fac=(-(LaRatio/SolidAngle)*(zta_o+1/2));
    
    MW_b1 = zeros(1,Nb); % Temp Variable; MultiThread
    MW_b2 = zeros(1,Nb); % Temp Variable; MultiThread
    
    for num_b = 1:Nb
        IntRange = [Node_eta_b(num_b)-detN_b Node_eta_b(num_b)];
        Int_Fn=@(eta,theta) HKI_Sub_SFn(1,eta,IntRange).*Int_W(eta,theta); % Define Integral Fn
        MW_b1(num_b) = HKI_Sub_Integral(Int_Fn,IntRange,IntError_b)*Fac;
        
        IntRange =  [Node_eta_b(num_b) Node_eta_b(num_b)+detN_b];
        Int_Fn=@(eta,theta) HKI_Sub_SFn(2,eta,IntRange).*Int_W(eta,theta); % Define Integral Fn
        MW_b2(num_b) = HKI_Sub_Integral(Int_Fn,IntRange,IntError_b)*Fac;
    end
    MW_b(num,:) = MW_b1 + MW_b2; % MD_b matrix construction; MultiThread
    
    %% Integral Ring, Surface 2, Global W
    R=@(zta,theta) HKI_Sub_R(eta_o,zta_o,1,theta,zta,LaRatio); % Define R
    Int_W=@(zta,theta) HKI_Sub_IntFn_W(2,R,ka,eta_o,0,theta,zta); % Define Integral Equation
    Fac=((LaRatio/SolidAngle));
    
    MW_r1 = zeros(1,Nr); % Temp Variable; MultiThread
    MW_r2 = zeros(1,Nr); % Temp Variable; MultiThread
    
    for num_r = 1:Nr
        IntRange = [Node_zta_r(num_r)-detN_r Node_zta_r(num_r)];
        Int_Fn=@(zta,theta) HKI_Sub_SFn(1,zta+0.5,IntRange+0.5).*Int_W(zta,theta); % Define Integral Fn
        MW_r1(num_r) = HKI_Sub_Integral(Int_Fn,IntRange,IntError_r)*Fac;
        
        IntRange =  [Node_zta_r(num_r) Node_zta_r(num_r)+detN_r];
        Int_Fn=@(zta,theta) HKI_Sub_SFn(2,zta+0.5,IntRange+0.5).*Int_W(zta,theta); % Define Integral Fn
        MW_r2(num_r) = HKI_Sub_Integral(Int_Fn,IntRange,IntError_r)*Fac;
    end
    MW_r(num,:) = MW_r1 + MW_r2; % MD_b matrix construction; MultiThread
    
    %% Integral Top, Surface 3, Global W
    R=@(eta,theta) HKI_Sub_R(eta_o,zta_o,eta,theta,0.5,LaRatio); % Define R
    Int_W=@(eta,theta) HKI_Sub_IntFn_W(3,R,ka,0,eta,theta,0); % Define Integral Equation
    Fac=((LaRatio/SolidAngle)*(zta_o-1/2));
    
    MW_t1 = zeros(1,Nt); % Temp Variable; MultiThread
    MW_t2 = zeros(1,Nt); % Temp Variable; MultiThread
    
    for num_t = 1:Nt
        IntRange = [Node_eta_t(num_t)-detN_t Node_eta_t(num_t)];
        Int_Fn=@(eta,theta) HKI_Sub_SFn(1,eta,IntRange).*Int_W(eta,theta); % Define Integral Fn
        MW_t1(num_t) = HKI_Sub_Integral(Int_Fn,IntRange,IntError_t)*Fac;
        
        IntRange =  [Node_eta_t(num_t) Node_eta_t(num_t)+detN_t];
        Int_Fn=@(eta,theta) HKI_Sub_SFn(2,eta,IntRange).*Int_W(eta,theta); % Define Integral Fn
        MW_t2(num_t) = HKI_Sub_Integral(Int_Fn,IntRange,IntError_t)*Fac;
    end
    MW_t(num,:) = MW_t1 + MW_t2; % MD_b matrix construction; MultiThread
end

%% Matrix Calculation
GD_ob(:,1:Nb) = MD_b;
GD_ob(:,Nb+1:Nb+Nr) = MD_r;
GD_ob(:,Nb+Nr+1:Nb+Nr+Nt) = MD_t;

GW_ob(:,1:length(MW_b(1,:))-1) = MW_b(:,1:length(MW_b(1,:))-1);
GW_ob(:,length(MW_b(1,:))) = MW_b(:,end) + MW_r(:,1);
GW_ob(:,length(MW_b(1,:))+1:(length(MW_b(1,:))+length(MW_r(1,:)))-2) = MW_r(:,2:1:length(MW_r(1,:))-1);
GW_ob(:,(length(MW_b(1,:))+length(MW_r(1,:)))-1) = MW_r(:,end) + MW_t(:,1);
GW_ob(:,(length(MW_b(1,:))+length(MW_r(1,:))):(length(MW_b(1,:))+length(MW_r(1,:))+length(MW_t(1,:)))-2) = MW_t(:,2:end);
end