function [GD, GW, p0] = HKI_Sub_SurPres_Mat(freq,u0,c0,rho0,radius_a,length_l,NN)

%% Geometry & Boundary Condition Parameter

Num = 2;
gap_g = 0.08;

total_length = Num*length_l + (Num-1)*gap_g;
totalLaRatio = total_length/radius_a;
L_totalL_Ratio = length_l/total_length;
g_totalL_Ratio = gap_g/total_length;

%% Define, number of node
Nb = NN(1); Nr = NN(2); Ng = NN(3); Nt = NN(end);

for i = 2 : 2*Num
    if mod(i,2) == 0
        NN(i) = Nr;
    elseif mod(i,2) == 1
        NN(i) = Ng;  
    end           
end
NN(2*Num+1) = Nt;

TotNumPre = sum(NN) - 2*Num; % Total number of Pressure Vector
TotNumVel = sum(NN); % Total number of Velocity Vector

detN_b = 1/(Nb-1); detN_r = L_totalL_Ratio/(Nr-1); detN_g = g_totalL_Ratio/(Ng-1); detN_t = 1/(Nt-1); % gap of nodes

% 적분 구간 설정시 특이점 제외를 위한 Error 정도 설정, 0.1%
IntError_b = detN_b*0.002; IntError_r = detN_r*0.001; IntError_g = detN_g*0.001; IntError_t = detN_t*0.002; 

%% Construct Node
% Node_eta_b = linspace(0,1,Nb)'; Node_eta_t = linspace(1,0,Nt)'; Node_zta_r = linspace(-0.5,0.5,Nr)';
Node_eta_b = (0:detN_b:1)'; Node_zta_r1 = (-0.5:detN_r:-0.5+L_totalL_Ratio)'; Node_zta_g = (-0.5+L_totalL_Ratio:detN_g:1/2*g_totalL_Ratio)'; Node_zta_r2 = (1/2*g_totalL_Ratio:detN_r:0.5)';Node_eta_t = (1:-detN_t:0)';

%% Construct node position
etaNode(1:NN(1)) = Node_eta_b;
ztaNode(1:NN(1)) = Node_zta_r1(1);

etaNode(NN(1)+1:sum(NN(1:2))-1) = Node_eta_b(end);
ztaNode(NN(1)+1:sum(NN(1:2))-1) = Node_zta_r1(2:end,1);

etaNode(sum(NN(1:2)):sum(NN(1:3))-2) = Node_eta_b(end);
ztaNode(sum(NN(1:2)):sum(NN(1:3))-2) = Node_zta_g(2:end,1);

etaNode(sum(NN(1:3))-1:sum(NN(1:4))-3) = Node_eta_b(end);
ztaNode(sum(NN(1:3))-1:sum(NN(1:4))-3) = Node_zta_r2(2:end,1);

etaNode(sum(NN(1:4))-2:sum(NN(1:5))-4) = Node_eta_t(2:end);
ztaNode(sum(NN(1:4))-2:sum(NN(1:5))-4) = Node_zta_r2(end);

%% Medium Parameter
% Water
c_water = c0; rho0_water = rho0;

ka = (2*pi*freq/c_water)*radius_a;
p0 = rho0_water*c_water*u0;

%% Cal Global D & W Matrix
% GD(PressureVectorNumber, VelocityVectorNumber)
GD = zeros(TotNumPre,TotNumVel);
GW = zeros(TotNumPre,TotNumPre);

MD_b = zeros(TotNumPre,Nb);
MD_ring1 = zeros(TotNumPre,Nr);
MD_gap = zeros(TotNumPre,Ng);
MD_ring2 = zeros(TotNumPre,Nr);
MD_t = zeros(TotNumPre,Nt);

MW_b = zeros(TotNumPre,Nb);
MW_ring1 = zeros(TotNumPre,Nr);
MW_gap = zeros(TotNumPre,Ng);
MW_ring2 = zeros(TotNumPre,Nr);
MW_t = zeros(TotNumPre,Nt);

parfor num = 1:TotNumPre
% for num = 1:TotNumPre
    % Define eto_o, zta_o
    eta_o = etaNode(num);
    zta_o = ztaNode(num);
    
    % Define Solid Angle
    if (num == Nb) || (num == Nb+Nr-1)
        SolidAngle = 3*pi;
    else
        SolidAngle = 2*pi;
    end
    
    %% Integral Bottom, Surface 1, Global D
    R=@(eta,theta) HKI_Sub_R(eta_o,zta_o,eta,theta,-0.5,totalLaRatio); % Define R using Sub_R(eta_o,zta_o,eta,theta,zta,LaRatio)
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
    
    %% Integral Ring 1, Surface 2, Global D
    R=@(zta,theta) HKI_Sub_R(eta_o,zta_o,1,theta,zta,totalLaRatio); % Define R
    Int_D=@(zta,theta) HKI_Sub_IntFn_D(2,R,ka,0,theta,zta); % Define Integral Equation
    Fac=((1i*ka)/SolidAngle)*totalLaRatio;
    
    MD_r1 = zeros(1,Nr); % Temp Variable; MultiThread
    MD_r2 = zeros(1,Nr); % Temp Variable; MultiThread
    
    for num_r = 1:Nr
        IntRange = [Node_zta_r1(num_r)-detN_r Node_zta_r1(num_r)];
        Int_Fn=@(zta,theta) HKI_Sub_SFn(1,zta+0.5,IntRange+0.5).*Int_D(zta,theta); % Define Integral Fn
        MD_r1(num_r) = HKI_Sub_Integral(Int_Fn,IntRange,IntError_r)*Fac;
             
        IntRange =  [Node_zta_r1(num_r) Node_zta_r1(num_r)+detN_r];
        Int_Fn=@(zta,theta) HKI_Sub_SFn(2,zta+0.5,IntRange+0.5).*Int_D(zta,theta); % Define Integral Fn
        MD_r2(num_r) = HKI_Sub_Integral(Int_Fn,IntRange,IntError_r)*Fac;
    
    end
    MD_ring1(num,:) = MD_r1 + MD_r2; % MD_b matrix construction; MultiThread

        %% Integral Gap, Surface 3, Global D
    R=@(zta,theta) HKI_Sub_R(eta_o,zta_o,1,theta,zta,totalLaRatio); % Define R
    Int_D=@(zta,theta) HKI_Sub_IntFn_D(2,R,ka,0,theta,zta); % Define Integral Equation
    Fac=((1i*ka)/SolidAngle)*totalLaRatio;
    
    MD_r1 = zeros(1,Ng); % Temp Variable; MultiThread
    MD_r2 = zeros(1,Ng); % Temp Variable; MultiThread
    
    for num_r = 1:Ng
        IntRange = [Node_zta_g(num_r)-detN_r Node_zta_g(num_r)];
        Int_Fn=@(zta,theta) HKI_Sub_SFn(1,zta+0.5,IntRange+0.5).*Int_D(zta,theta); % Define Integral Fn
        MD_r1(num_r) = HKI_Sub_Integral(Int_Fn,IntRange,IntError_r)*Fac;
             
        IntRange =  [Node_zta_g(num_r) Node_zta_g(num_r)+detN_r];
        Int_Fn=@(zta,theta) HKI_Sub_SFn(2,zta+0.5,IntRange+0.5).*Int_D(zta,theta); % Define Integral Fn
        MD_r2(num_r) = HKI_Sub_Integral(Int_Fn,IntRange,IntError_r)*Fac;
    
    end
    MD_gap(num,:) = MD_r1 + MD_r2; % MD_b matrix construction; MultiThread

        %% Integral Ring 2, Surface 4, Global D
    R=@(zta,theta) HKI_Sub_R(eta_o,zta_o,1,theta,zta,totalLaRatio); % Define R
    Int_D=@(zta,theta) HKI_Sub_IntFn_D(2,R,ka,0,theta,zta); % Define Integral Equation
    Fac=((1i*ka)/SolidAngle)*totalLaRatio;
    
    MD_r1 = zeros(1,Nr); % Temp Variable; MultiThread
    MD_r2 = zeros(1,Nr); % Temp Variable; MultiThread
    
    for num_r = 1:Nr
        IntRange = [Node_zta_r2(num_r)-detN_r Node_zta_r2(num_r)];
        Int_Fn=@(zta,theta) HKI_Sub_SFn(1,zta+0.5,IntRange+0.5).*Int_D(zta,theta); % Define Integral Fn
        MD_r1(num_r) = HKI_Sub_Integral(Int_Fn,IntRange,IntError_r)*Fac;
             
        IntRange =  [Node_zta_r2(num_r) Node_zta_r2(num_r)+detN_r];
        Int_Fn=@(zta,theta) HKI_Sub_SFn(2,zta+0.5,IntRange+0.5).*Int_D(zta,theta); % Define Integral Fn
        MD_r2(num_r) = HKI_Sub_Integral(Int_Fn,IntRange,IntError_r)*Fac;
    
    end
    MD_ring2(num,:) = MD_r1 + MD_r2; % MD_b matrix construction; MultiThread
    
    %% Integral Top, Surface 5, Global D
    R=@(eta,theta) HKI_Sub_R(eta_o,zta_o,eta,theta,0.5,totalLaRatio); % Define R
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
    R=@(eta,theta) HKI_Sub_R(eta_o,zta_o,eta,theta,-0.5,totalLaRatio); % Define R
    Int_W=@(eta,theta) HKI_Sub_IntFn_W(1,R,ka,0,eta,theta,0); % Define Integral Equation
    Fac=(-(totalLaRatio/SolidAngle)*(zta_o+1/2));
    
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
    
    %% Integral Ring 1, Surface 2, Global W
    R=@(zta,theta) HKI_Sub_R(eta_o,zta_o,1,theta,zta,totalLaRatio); % Define R
    Int_W=@(zta,theta) HKI_Sub_IntFn_W(2,R,ka,eta_o,0,theta,zta); % Define Integral Equation
    Fac=((totalLaRatio/SolidAngle));
    
    MW_r1 = zeros(1,Nr); % Temp Variable; MultiThread
    MW_r2 = zeros(1,Nr); % Temp Variable; MultiThread
    
    for num_r = 1:Nr
        IntRange = [Node_zta_r1(num_r)-detN_r Node_zta_r1(num_r)];
        Int_Fn=@(zta,theta) HKI_Sub_SFn(1,zta+0.5,IntRange+0.5).*Int_W(zta,theta); % Define Integral Fn
        MW_r1(num_r) = HKI_Sub_Integral(Int_Fn,IntRange,IntError_r)*Fac;
                
        IntRange =  [Node_zta_r1(num_r) Node_zta_r1(num_r)+detN_r];
        Int_Fn=@(zta,theta) HKI_Sub_SFn(2,zta+0.5,IntRange+0.5).*Int_W(zta,theta); % Define Integral Fn
        MW_r2(num_r) = HKI_Sub_Integral(Int_Fn,IntRange,IntError_r)*Fac;
        
    end
    MW_ring1(num,:) = MW_r1 + MW_r2; % MD_b matrix construction; MultiThread

        %% Integral Gap, Surface 2, Global W
    R=@(zta,theta) HKI_Sub_R(eta_o,zta_o,1,theta,zta,totalLaRatio); % Define R
    Int_W=@(zta,theta) HKI_Sub_IntFn_W(2,R,ka,eta_o,0,theta,zta); % Define Integral Equation
    Fac=((totalLaRatio/SolidAngle));
    
    MW_r1 = zeros(1,Ng); % Temp Variable; MultiThread
    MW_r2 = zeros(1,Ng); % Temp Variable; MultiThread
    
    for num_r = 1:Ng
        IntRange = [Node_zta_g(num_r)-detN_r Node_zta_g(num_r)];
        Int_Fn=@(zta,theta) HKI_Sub_SFn(1,zta+0.5,IntRange+0.5).*Int_W(zta,theta); % Define Integral Fn
        MW_r1(num_r) = HKI_Sub_Integral(Int_Fn,IntRange,IntError_r)*Fac;
                
        IntRange =  [Node_zta_g(num_r) Node_zta_g(num_r)+detN_r];
        Int_Fn=@(zta,theta) HKI_Sub_SFn(2,zta+0.5,IntRange+0.5).*Int_W(zta,theta); % Define Integral Fn
        MW_r2(num_r) = HKI_Sub_Integral(Int_Fn,IntRange,IntError_r)*Fac;
        
    end
    MW_gap(num,:) = MW_r1 + MW_r2; % MD_b matrix construction; MultiThread

        %% Integral Ring, Surface 2, Global W
    R=@(zta,theta) HKI_Sub_R(eta_o,zta_o,1,theta,zta,totalLaRatio); % Define R
    Int_W=@(zta,theta) HKI_Sub_IntFn_W(2,R,ka,eta_o,0,theta,zta); % Define Integral Equation
    Fac=((totalLaRatio/SolidAngle));
    
    MW_r1 = zeros(1,Nr); % Temp Variable; MultiThread
    MW_r2 = zeros(1,Nr); % Temp Variable; MultiThread
    
    for num_r = 1:Nr
        IntRange = [Node_zta_r2(num_r)-detN_r Node_zta_r2(num_r)];
        Int_Fn=@(zta,theta) HKI_Sub_SFn(1,zta+0.5,IntRange+0.5).*Int_W(zta,theta); % Define Integral Fn
        MW_r1(num_r) = HKI_Sub_Integral(Int_Fn,IntRange,IntError_r)*Fac;
                
        IntRange =  [Node_zta_r2(num_r) Node_zta_r2(num_r)+detN_r];
        Int_Fn=@(zta,theta) HKI_Sub_SFn(2,zta+0.5,IntRange+0.5).*Int_W(zta,theta); % Define Integral Fn
        MW_r2(num_r) = HKI_Sub_Integral(Int_Fn,IntRange,IntError_r)*Fac;
        
    end
    MW_ring2(num,:) = MW_r1 + MW_r2; % MD_b matrix construction; MultiThread
    %% Integral Top, Surface 3, Global W
    R=@(eta,theta) HKI_Sub_R(eta_o,zta_o,eta,theta,0.5,totalLaRatio); % Define R
    Int_W=@(eta,theta) HKI_Sub_IntFn_W(3,R,ka,0,eta,theta,0); % Define Integral Equation
    Fac=((totalLaRatio/SolidAngle)*(zta_o-1/2));
    
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
GD(:,1:NN(1)) = MD_b;
GD(:,NN(1)+1:sum(NN(1:2))) = MD_ring1;
GD(:,sum(NN(1:2))+1:sum(NN(1:3))) = MD_gap;
GD(:,sum(NN(1:3))+1:sum(NN(1:4))) = MD_ring2;
GD(:,sum(NN(1:4))+1:sum(NN(1:5))) = MD_t;

GW(:,1:NN(1)-1) = MW_b(:,1:NN(1)-1);
GW(:,NN(1)) = MW_b(:,end) + MW_ring1(:,1);
GW(:,NN(1)+1:sum(NN(1:2))-2) = MW_ring1(:,2:NN(2)-1);

GW(:,sum(NN(1:2))-1) = MW_ring1(:,end) + MW_gap(:,1);
GW(:,sum(NN(1:2)):sum(NN(1:3))-3) = MW_gap(:,2:NN(3)-1);

GW(:,sum(NN(1:3))-2) = MW_gap(:,end) + MW_ring2(:,1);
GW(:,sum(NN(1:3))-1:sum(NN(1:4))-4) = MW_ring2(:,2:NN(4)-1);

GW(:,sum(NN(1:4))-3) = MW_ring2(:,end) + MW_t(:,1);
GW(:,sum(NN(1:4))-2:sum(NN(1:5))-4) = MW_t(:,2:end);

end



























