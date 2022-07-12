function [z_r]=RadiationImp_ParameterStudy_4Array_Subrutine(rhowater,sp,ka,a,L,RL,g_side,g_center,zlength,RadiationTable,TableParameter)
% Calculate z_sR [Page 48]
% {1,1} = t_t, {1.2} = t_r1, {1,3} = t_g1, {1,4} = t_r2, {1,5} = t_g2,
% {1,6} = t_r3, {1,7} = t_g3, {1,8} = t_r4, {1,9} = t_b,
% {1,10} = r1_r1, {1.11} = r1_g1, {1,12} = r1_r2, {1,13} = r1_g2, 
% {1,14} = r1_r3, {1,15} = r1_g3, {1,16} = r1_r4,
% {1,17} = g1_g1, {1,18} = g1_r2, {1,19} = g1_g2, {1,20} = g1_r3,
% {1,21} = g1_g3, 
% {1,22} = r2_r2, {1,23} = r2_g2, {1,24} = r2_r3, {1,25} = g2_g2,
% {1,26} = r3_r3, {1,27} = g3_g3, {1,28} = r4_r4.

%% Sort Impedance data using Table
Meshka = TableParameter{1,1};
Meshga = TableParameter{1,2};
TablegaRatio = TableParameter{1,3};

% 기존 데이터로 만들기 위한 dimensionless factor
MultiFactor_L = rhowater*sp*2*pi*a*RL;
MultiFactor_t = rhowater*sp*pi*a^2;
MultiFactor_g_side = rhowater*sp*2*pi*a*g_side;
MultiFactor_g_center = rhowater*sp*2*pi*a*g_center;

ZttRealData = real(RadiationTable{1,1});
ZttImagData = imag(RadiationTable{1,1});

Ztr1RealData = real(RadiationTable{1,2});
Ztr1ImagData = imag(RadiationTable{1,2});

Ztg1RealData = real(RadiationTable{1,3});
Ztg1ImagData = imag(RadiationTable{1,3});

Ztr2RealData = real(RadiationTable{1,4});
Ztr2ImagData = imag(RadiationTable{1,4});

Ztg2RealData = real(RadiationTable{1,5});
Ztg2ImagData = imag(RadiationTable{1,5});

Ztr3RealData = real(RadiationTable{1,6});
Ztr3ImagData = imag(RadiationTable{1,6});

Ztg3RealData = real(RadiationTable{1,7});
Ztg3ImagData = imag(RadiationTable{1,7});

Ztr4RealData = real(RadiationTable{1,8});
Ztr4ImagData = imag(RadiationTable{1,8});

ZtbRealData = real(RadiationTable{1,9});
ZtbImagData = imag(RadiationTable{1,9});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Zr1r1RealData = real(RadiationTable{1,10});
Zr1r1ImagData = imag(RadiationTable{1,10});

Zr1g1RealData = real(RadiationTable{1,11});
Zr1g1ImagData = imag(RadiationTable{1,11});

Zr1r2RealData = real(RadiationTable{1,12});
Zr1r2ImagData = imag(RadiationTable{1,12});

Zr1g2RealData = real(RadiationTable{1,13});
Zr1g2ImagData = imag(RadiationTable{1,13});

Zr1r3RealData = real(RadiationTable{1,14});
Zr1r3ImagData = imag(RadiationTable{1,14});

Zr1g3RealData = real(RadiationTable{1,15});
Zr1g3ImagData = imag(RadiationTable{1,15});

Zr1r4RealData = real(RadiationTable{1,16});
Zr1r4ImagData = imag(RadiationTable{1,16});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Zg1g1RealData = real(RadiationTable{1,17});
Zg1g1ImagData = imag(RadiationTable{1,17});

Zg1r2RealData = real(RadiationTable{1,18});
Zg1r2ImagData = imag(RadiationTable{1,18});

Zg1g2RealData = real(RadiationTable{1,19});
Zg1g2ImagData = imag(RadiationTable{1,19});

Zg1r3RealData = real(RadiationTable{1,20});
Zg1r3ImagData = imag(RadiationTable{1,20});

Zg1g3RealData = real(RadiationTable{1,21});
Zg1g3ImagData = imag(RadiationTable{1,21});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Zr2r2RealData = real(RadiationTable{1,22});
Zr2r2ImagData = imag(RadiationTable{1,22});

Zr2g2RealData = real(RadiationTable{1,23});
Zr2g2ImagData = imag(RadiationTable{1,23});

Zr2r3RealData = real(RadiationTable{1,24});
Zr2r3ImagData = imag(RadiationTable{1,24});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Zg2g2RealData = real(RadiationTable{1,25});
Zg2g2ImagData = imag(RadiationTable{1,25});

Zr3r3RealData = real(RadiationTable{1,26});
Zr3r3ImagData = imag(RadiationTable{1,26});

Zg3g3RealData = real(RadiationTable{1,27});
Zg3g3ImagData = imag(RadiationTable{1,27});

Zr4r4RealData = real(RadiationTable{1,28});
Zr4r4ImagData = imag(RadiationTable{1,28});

%% Find Radiation Impedance Data using Table for FFR ring

Fac = g_side/a; % gaRatio

RealZtt=interp2(Meshka,Meshga,ZttRealData,ka,Fac,'spline');
ImagZtt=interp2(Meshka,Meshga,ZttImagData,ka,Fac,'spline');

RealZtr1=interp2(Meshka,Meshga,Ztr1RealData,ka,Fac,'spline');
ImagZtr1=interp2(Meshka,Meshga,Ztr1ImagData,ka,Fac,'spline');

RealZtg1=interp2(Meshka,Meshga,Ztg1RealData,ka,Fac,'spline');
ImagZtg1=interp2(Meshka,Meshga,Ztg1ImagData,ka,Fac,'spline');

RealZtr2=interp2(Meshka,Meshga,Ztr2RealData,ka,Fac,'spline');
ImagZtr2=interp2(Meshka,Meshga,Ztr2ImagData,ka,Fac,'spline');

RealZtg2=interp2(Meshka,Meshga,Ztg2RealData,ka,Fac,'spline');
ImagZtg2=interp2(Meshka,Meshga,Ztg2ImagData,ka,Fac,'spline');

RealZtr3=interp2(Meshka,Meshga,Ztr3RealData,ka,Fac,'spline');
ImagZtr3=interp2(Meshka,Meshga,Ztr3ImagData,ka,Fac,'spline');

RealZtg3=interp2(Meshka,Meshga,Ztg3RealData,ka,Fac,'spline');
ImagZtg3=interp2(Meshka,Meshga,Ztg3ImagData,ka,Fac,'spline');

RealZtr4=interp2(Meshka,Meshga,Ztr4RealData,ka,Fac,'spline');
ImagZtr4=interp2(Meshka,Meshga,Ztr4ImagData,ka,Fac,'spline');

RealZtb=interp2(Meshka,Meshga,ZtbRealData,ka,Fac,'spline');
ImagZtb=interp2(Meshka,Meshga,ZtbImagData,ka,Fac,'spline');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RealZr1r1=interp2(Meshka,Meshga,Zr1r1RealData,ka,Fac,'spline');
ImagZr1r1=interp2(Meshka,Meshga,Zr1r1ImagData,ka,Fac,'spline');

RealZr1g1=interp2(Meshka,Meshga,Zr1g1RealData,ka,Fac,'spline');
ImagZr1g1=interp2(Meshka,Meshga,Zr1g1ImagData,ka,Fac,'spline');

RealZr1r2=interp2(Meshka,Meshga,Zr1r2RealData,ka,Fac,'spline');
ImagZr1r2=interp2(Meshka,Meshga,Zr1r2ImagData,ka,Fac,'spline');

RealZr1g2=interp2(Meshka,Meshga,Zr1g2RealData,ka,Fac,'spline');
ImagZr1g2=interp2(Meshka,Meshga,Zr1g2ImagData,ka,Fac,'spline');

RealZr1r3=interp2(Meshka,Meshga,Zr1r3RealData,ka,Fac,'spline');
ImagZr1r3=interp2(Meshka,Meshga,Zr1r3ImagData,ka,Fac,'spline');

RealZr1g3=interp2(Meshka,Meshga,Zr1g3RealData,ka,Fac,'spline');
ImagZr1g3=interp2(Meshka,Meshga,Zr1g3ImagData,ka,Fac,'spline');

RealZr1r4=interp2(Meshka,Meshga,Zr1r4RealData,ka,Fac,'spline');
ImagZr1r4=interp2(Meshka,Meshga,Zr1r4ImagData,ka,Fac,'spline');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RealZg1g1=interp2(Meshka,Meshga,Zg1g1RealData,ka,Fac,'spline');
ImagZg1g1=interp2(Meshka,Meshga,Zg1g1ImagData,ka,Fac,'spline');

RealZg1r2=interp2(Meshka,Meshga,Zg1r2RealData,ka,Fac,'spline');
ImagZg1r2=interp2(Meshka,Meshga,Zg1r2ImagData,ka,Fac,'spline');

RealZg1g2=interp2(Meshka,Meshga,Zg1g2RealData,ka,Fac,'spline');
ImagZg1g2=interp2(Meshka,Meshga,Zg1g2ImagData,ka,Fac,'spline');

RealZg1r3=interp2(Meshka,Meshga,Zg1r3RealData,ka,Fac,'spline');
ImagZg1r3=interp2(Meshka,Meshga,Zg1r3ImagData,ka,Fac,'spline');

RealZg1g3=interp2(Meshka,Meshga,Zg1g3RealData,ka,Fac,'spline');
ImagZg1g3=interp2(Meshka,Meshga,Zg1g3ImagData,ka,Fac,'spline');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RealZr2r2=interp2(Meshka,Meshga,Zr2r2RealData,ka,Fac,'spline');
ImagZr2r2=interp2(Meshka,Meshga,Zr2r2ImagData,ka,Fac,'spline');

RealZr2g2=interp2(Meshka,Meshga,Zr2g2RealData,ka,Fac,'spline');
ImagZr2g2=interp2(Meshka,Meshga,Zr2g2ImagData,ka,Fac,'spline');

RealZr2r3=interp2(Meshka,Meshga,Zr2r3RealData,ka,Fac,'spline');
ImagZr2r3=interp2(Meshka,Meshga,Zr2r3ImagData,ka,Fac,'spline');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RealZg2g2=interp2(Meshka,Meshga,Zg2g2RealData,ka,Fac,'spline');
ImagZg2g2=interp2(Meshka,Meshga,Zg2g2ImagData,ka,Fac,'spline');

RealZr3r3=interp2(Meshka,Meshga,Zr3r3RealData,ka,Fac,'spline');
ImagZr3r3=interp2(Meshka,Meshga,Zr3r3ImagData,ka,Fac,'spline');

RealZg3g3=interp2(Meshka,Meshga,Zg3g3RealData,ka,Fac,'spline');
ImagZg3g3=interp2(Meshka,Meshga,Zg3g3ImagData,ka,Fac,'spline');

RealZr4r4=interp2(Meshka,Meshga,Zr4r4RealData,ka,Fac,'spline');
ImagZr4r4=interp2(Meshka,Meshga,Zr4r4ImagData,ka,Fac,'spline');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ztt = (RealZtt + ImagZtt*1i) * MultiFactor_t;
Ztr1 = (RealZtr1 + ImagZtr1*1i) * MultiFactor_t;
Ztg1 = (RealZtg1 + ImagZtg1*1i) * MultiFactor_t;
Ztr2 = (RealZtr2 + ImagZtr2*1i) * MultiFactor_t;
Ztg2 = (RealZtg2 + ImagZtg2*1i) * MultiFactor_t;
Ztr3 = (RealZtr3 + ImagZtr3*1i) * MultiFactor_t;
Ztg3 = (RealZtg3 + ImagZtg3*1i) * MultiFactor_t;
Ztr4 = (RealZtr4 + ImagZtr4*1i) * MultiFactor_t;
Ztb = (RealZtb + ImagZtb*1i) * MultiFactor_t;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Zr1r1 = (RealZr1r1 + ImagZr1r1*1i) * MultiFactor_L;
Zr1g1 = (RealZr1g1 + ImagZr1g1*1i) * MultiFactor_L;
Zr1r2 = (RealZr1r2 + ImagZr1r2*1i) * MultiFactor_L;
Zr1g2 = (RealZr1g2 + ImagZr1g2*1i) * MultiFactor_L;
Zr1r3 = (RealZr1r3 + ImagZr1r3*1i) * MultiFactor_L;
Zr1g3 = (RealZr1g3 + ImagZr1g3*1i) * MultiFactor_L;
Zr1r4 = (RealZr1r4 + ImagZr1r4*1i) * MultiFactor_L;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Zg1g1 = (RealZg1g1 + ImagZg1g1*1i) * MultiFactor_g_side;
Zg1r2 = (RealZg1r2 + ImagZg1r2*1i) * MultiFactor_g_side;
Zg1g2 = (RealZg1g2 + ImagZg1g2*1i) * MultiFactor_g_side;
Zg1r3 = (RealZg1r3 + ImagZg1r3*1i) * MultiFactor_g_side;
Zg1g3 = (RealZg1g3 + ImagZg1g3*1i) * MultiFactor_g_side;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Zr2r2 = (RealZr2r2 + ImagZr2r2*1i) * MultiFactor_L;
Zr2g2 = (RealZr2g2 + ImagZr2g2*1i) * MultiFactor_L;
Zr2r3 = (RealZr2r3 + ImagZr2r3*1i) * MultiFactor_L;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Zg2g2 = (RealZg2g2 + ImagZg2g2*1i) * MultiFactor_g_center;
Zr3r3 = (RealZr3r3 + ImagZr3r3*1i) * MultiFactor_L;
Zg3g3 = (RealZg3g3 + ImagZg3g3*1i) * MultiFactor_g_side;
Zr4r4 = (RealZr4r4 + ImagZr4r4*1i) * MultiFactor_L;

% Import Data from Excel File
% sheet numbers -> %

%% Construce Radiation Impednace Matrix
z_r=cell(1,length(ka)); % radiation impedance matrix

for i=1:length(ka);
    z_r{1,i}=zeros(zlength,zlength);

    %% 4Array

    % from Top
    z_r{1,i}(1,1)=Ztt(i)/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(1,2)=Ztr1(i)/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(1,4)=Ztg1(i)/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(1,6)=Ztr2(i)/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(1,8)=Ztg2(i)/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(1,10)=Ztr3(i)/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(1,12)=Ztg3(i)/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(1,14)=Ztr4(i)/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(1,15)=-Ztb(i)/(rhowater*sp*2*pi*a*L);
    
    % from 1st ring (Top ring)
    z_r{1,i}(2,1)=z_r{1,i}(1,2);
    z_r{1,i}(2,2)=Zr1r1(i)/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(2,4)=Zr1g1(i)/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(2,6)=Zr1r2(i)/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(2,8)=Zr1g2(i)/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(2,10)=Zr1r3(i)/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(2,12)=Zr1g3(i)/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(2,14)=Zr1r4(i)/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(2,15)=-z_r{1,i}(1,14);
    
    % from 1st gab
    z_r{1,i}(4,1)=z_r{1,i}(1,4);
    z_r{1,i}(4,2)=z_r{1,i}(2,4);
    z_r{1,i}(4,4)=Zg1g1(i)/(rhowater*sp*2*pi*a*L); 
    z_r{1,i}(4,6)=Zg1r2(i)/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(4,8)=Zg1g2(i)/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(4,10)=Zg1r3(i)/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(4,12)=Zg1g3(i)/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(4,14)=z_r{1,i}(2,12);
    z_r{1,i}(4,15)=-z_r{1,i}(1,12);
    
    % from 2nd ring (Mid 1 ring)
    z_r{1,i}(6,1)=z_r{1,i}(1,6);
    z_r{1,i}(6,2)=z_r{1,i}(2,6);
    z_r{1,i}(6,4)=z_r{1,i}(4,6); 
    z_r{1,i}(6,6)=Zr2r2(i)/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(6,8)=Zr2g2(i)/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(6,10)=Zr2r3(i)/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(6,12)=z_r{1,i}(4,10);
    z_r{1,i}(6,14)=z_r{1,i}(2,10);
    z_r{1,i}(6,15)=-z_r{1,i}(1,10);
    
    % from 2nd gab
    z_r{1,i}(8,1)=z_r{1,i}(1,8);
    z_r{1,i}(8,2)=z_r{1,i}(2,8);
    z_r{1,i}(8,4)=z_r{1,i}(4,8); 
    z_r{1,i}(8,6)=z_r{1,i}(6,8);
    z_r{1,i}(8,8)=Zg2g2(i)/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(8,10)=z_r{1,i}(6,8);
    z_r{1,i}(8,12)=z_r{1,i}(4,8);
    z_r{1,i}(8,14)=z_r{1,i}(2,8);
    z_r{1,i}(8,15)=-z_r{1,i}(1,8);
    
    % from 3nd ring (Mid 2 ring)
    z_r{1,i}(10,1)=z_r{1,i}(1,10);
    z_r{1,i}(10,2)=z_r{1,i}(2,10);
    z_r{1,i}(10,4)=z_r{1,i}(4,10); 
    z_r{1,i}(10,6)=z_r{1,i}(6,10);
    z_r{1,i}(10,8)=z_r{1,i}(8,10);
    z_r{1,i}(10,10)=Zr3r3(i)/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(10,12)=z_r{1,i}(4,6);
    z_r{1,i}(10,14)=z_r{1,i}(2,6);
    z_r{1,i}(10,15)=-z_r{1,i}(1,6);

    % from 3rd gab 
    z_r{1,i}(12,1)=z_r{1,i}(1,12);
    z_r{1,i}(12,2)=z_r{1,i}(2,12);
    z_r{1,i}(12,4)=z_r{1,i}(4,12); 
    z_r{1,i}(12,6)=z_r{1,i}(6,12);
    z_r{1,i}(12,8)=z_r{1,i}(8,12);
    z_r{1,i}(12,10)=z_r{1,i}(10,12);
    z_r{1,i}(12,12)=Zg3g3(i)/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(12,14)=z_r{1,i}(2,4);
    z_r{1,i}(12,15)=-z_r{1,i}(1,4);
    
    % from 4th ring
    z_r{1,i}(14,1)=z_r{1,i}(1,14);
    z_r{1,i}(14,2)=z_r{1,i}(2,14);
    z_r{1,i}(14,4)=z_r{1,i}(4,14); 
    z_r{1,i}(14,6)=z_r{1,i}(6,14);
    z_r{1,i}(14,8)=z_r{1,i}(8,14);
    z_r{1,i}(14,10)=z_r{1,i}(10,14);
    z_r{1,i}(14,12)=z_r{1,i}(12,14);
    z_r{1,i}(14,14)=Zr4r4(i)/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(14,15)=-z_r{1,i}(1,2);
    
    % from bottom
    z_r{1,i}(15,1)=z_r{1,i}(1,15);
    z_r{1,i}(15,2)=z_r{1,i}(2,15);
    z_r{1,i}(15,4)=z_r{1,i}(4,15); 
    z_r{1,i}(15,6)=z_r{1,i}(6,15);
    z_r{1,i}(15,8)=z_r{1,i}(8,15);
    z_r{1,i}(15,10)=z_r{1,i}(10,15);
    z_r{1,i}(15,12)=z_r{1,i}(12,15);
    z_r{1,i}(15,14)=z_r{1,i}(14,15);
    z_r{1,i}(15,15)=z_r{1,i}(1,1);
        
    %% For circuit model
%     % from 1st ring
%     z_r{1,i}(2,1)=0;
%     z_r{1,i}(2,2)=exceldata{1,2}{i,1}/(rhowater*sp*2*pi*a*L);
%     z_r{1,i}(2,3)=0;
% 
%     %from bottom, top
%     z_r{1,i}(1,1)=exceldata{1,3}{i,1}/(rhowater*sp*2*pi*a*L);
%     z_r{1,i}(1,2)=0;
%     z_r{1,i}(1,3)=0;
%     
%     z_r{1,i}(3,1)=0;
%     z_r{1,i}(3,2)=0;
%     z_r{1,i}(3,3)=z_r{1,i}(1,1);   

end
end