function [ExportCircuit,UMatrix_for_HKI,Time_Main_Sub_Single_FFR] = Main_Sub_Single_FFR(freq,c0,rho0,radius_a,length_l,z_RMatrix)
%% FFR Circuit Model Code
% Made by Seungwon Nam
% Modified by Kyounghun Been

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
rhowater=rho0; % Density of Water
sp=c0; % Sound Speed at Water

f = freq;
omega=2*pi*f; % Angular Frequency
k=@(w) w/sp; % Wave Number

num=1;
zlength=3+4*(num-1); % length of the impedance matrix of array number = num

a=radius_a; % Ring Average Radius
L=length_l; % Ring Height
t=0.04; % Ring Thickness

LaRatio = L/a; % L/a ratio, dimensionless parameter by Been
ka = ((2*pi.*f)./sp).*a; % ka, dimensionless parameter by Been
Dimfactor = rhowater*sp*2*pi*a*L;

n=500; % Sum Number (n>50 is enough)

%% Circuit Parameters of PZT
% Ignore Electrical Dissipation G0 and Mechanical Damping Rm

C0=(2*pi*a*L/t)*eS33; % Equivalent Clamped Capacitance
N=(2*pi*L)*(d31/sE11); % Equivalent Electromechanical turns Ratio
M=rhopzt*2*pi*a*L*t; % Mass of PZT
CE=(sE11*a)/(2*pi*t*L); % Equivalent Compliance

f0=1/(2*pi)*sqrt((N^2+C0/CE)/(M*C0)); % Resonance Frequency of Uncoupled PZT Ring
ka0 = ((2*pi*f0)/sp)*a; % ka of Uncoupled PZT Ring
kaRatio = ka/ka0;

%% Find the Lists of Zeros of J_0'(x)=0
% Checked accuracy until n=10000

J0dot=@(x) -besselj(1,x);
rootj0dot=zeros(1,length(f));

for i=1:n;
    rootj0dot(i)=fzero(J0dot,[i i+1]*pi); % Need to ignore first zero x=0
end

%% Circuit Parameters of Inner Fluid
% 검토 완료

[z_rr,z_rz,z_1,z_2,z_PRm]=FFR_Sub_CavImp(LaRatio,t,a,ka,rhopzt,rhowater,sp,sE11,n,rootj0dot);
Z_rr = z_rr.*Dimfactor;
Z_rz = z_rz.*Dimfactor;
Z_1 = z_1.*Dimfactor;
Z_2 = z_2.*Dimfactor;
Z_PRm = z_PRm.*Dimfactor;

%% Circuit Parameters about Radiation Impedance (For 4Array)

% [z_RMatrix, ~]=FFR_Sub_RadImp(VF_For_FFR,rhowater,sp,ka,a,L,zlength);

%% Calculate Total Admittance
clear Y;

% Construct Matrix and Vector
z_cMatrix = cell(1,length(ka)); % Cavity mode impedance matrix
z_PRMatrix = cell(1,length(ka)); % Piezoelectric ring matrix
for i=1:length(ka);
    z_cMatrix{1,i} = [z_1(i), -z_rz(i), -z_2(i); -z_rz(i), z_rr(i), z_rz(i); -z_2(i), z_rz(i), z_1(i)];
    z_PRMatrix{1,i} = [0 0 0; 0 z_PRm(i) 0; 0 0 0];
end
VectorT2 = [0;1;0];

Y = zeros(1,length(ka));
Admittance_Factor = N^2/(rhowater*sp*2*pi*a*L);
InversMatrix= cell(1,length(ka));
for i=1:length(ka);
    InversMatrix{1,i} =(inv(z_cMatrix{1,i} + z_PRMatrix{1,i} + z_RMatrix{i,1}));
    Y_temp = VectorT2' * InversMatrix{1,i} * VectorT2;
    Y(i) = (1i*(2*pi*f(i))*C0) + Admittance_Factor * Y_temp ;
end

ExportCircuit = [f', kaRatio', Y.'];

%% Temp
% uFactor = d31/(a*sE11*rhowater*sp^2);
uFactor = d31/(a*sE11);
uMatrix = cell(1,length(ka));
UMatrix_for_HKI = cell(1,length(ka));

for i=1:length(ka);
    u_temp = uFactor.*InversMatrix{1,i}*VectorT2;
    uMatrix{1,i} = u_temp;
    UMatrix_for_HKI{1,i} = u_temp;
end

Time_Main_Sub_Single_FFR = toc;

Hour = fix(Time_Main_Sub_Single_FFR/3600); Min = fix(rem(Time_Main_Sub_Single_FFR,3600)/60); Sec = round(rem(rem(Time_Main_Sub_Single_FFR,3600),60));
fprintf('Main_Sub_Single_FFR 계산 소요 시간은 %d시간 %d분 %d초 입니다.\n',Hour,Min,Sec)
end
