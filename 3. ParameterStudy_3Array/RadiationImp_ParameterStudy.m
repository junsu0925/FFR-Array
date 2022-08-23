function [z_r]=RadiationImp_ParameterStudy(rhowater,sp,ka,a,L,g,zlength,RadiationTable,TableParameter)
% Calculate z_sR [Page 48]
% {1,1} = r1_b, {1.2} = r1_r1, {1,3} = r1_g, {1,4} = r1_r2, {1,5} = r1_t
% {1,6} = g_b {1,7} = g_g1 {1,8} = b_b {1,9} = b_t

%% Sort Impedance data using Table
Meshka = TableParameter{1,1};
Meshga = TableParameter{1,2};
TablegaRatio = TableParameter{1,3};

% 기존 데이터로 만들기 위한 dimensionless factor
MultiFactor_L = rhowater*sp*2*pi*a*L;
MultiFactor_b = rhowater*sp*pi*a^2;
MultiFactor_g = rhowater*sp*2*pi*a*g;

Zr1bRealData = real(RadiationTable{1,1});
Zr1bImagData = imag(RadiationTable{1,1});

Zr1r1RealData = real(RadiationTable{1,2});
Zr1r1ImagData = imag(RadiationTable{1,2});

Zr1gRealData = real(RadiationTable{1,3});
Zr1gImagData = imag(RadiationTable{1,3});

Zr1r2RealData = real(RadiationTable{1,4});
Zr1r2ImagData = imag(RadiationTable{1,4});

Zr1tRealData = real(RadiationTable{1,5});
Zr1tImagData = imag(RadiationTable{1,5});

ZgbRealData = real(RadiationTable{1,6});
ZgbImagData = imag(RadiationTable{1,6});

ZggRealData = real(RadiationTable{1,7});
ZggImagData = imag(RadiationTable{1,7});

ZbbRealData = real(RadiationTable{1,8});
ZbbImagData = imag(RadiationTable{1,8});

ZbtRealData = real(RadiationTable{1,9});
ZbtImagData = imag(RadiationTable{1,9});

%% Find Radiation Impedance Data using Table for FFR ring

Fac = g/a; % gaRatio

RealZr1b=interp2(Meshka,Meshga,Zr1bRealData,ka,Fac,'spline');
ImagZr1b=interp2(Meshka,Meshga,Zr1bImagData,ka,Fac,'spline');

RealZr1r1=interp2(Meshka,Meshga,Zr1r1RealData,ka,Fac,'spline');
ImagZr1r1=interp2(Meshka,Meshga,Zr1r1ImagData,ka,Fac,'spline');

RealZr1g=interp2(Meshka,Meshga,Zr1gRealData,ka,Fac,'spline');
ImagZr1g=interp2(Meshka,Meshga,Zr1gImagData,ka,Fac,'spline');

RealZr1r2=interp2(Meshka,Meshga,Zr1r2RealData,ka,Fac,'spline');
ImagZr1r2=interp2(Meshka,Meshga,Zr1r2ImagData,ka,Fac,'spline');

RealZr1t=interp2(Meshka,Meshga,Zr1tRealData,ka,Fac,'spline');
ImagZr1t=interp2(Meshka,Meshga,Zr1tImagData,ka,Fac,'spline');

RealZgb=interp2(Meshka,Meshga,ZgbRealData,ka,Fac,'spline');
ImagZgb=interp2(Meshka,Meshga,ZgbImagData,ka,Fac,'spline');

RealZgg=interp2(Meshka,Meshga,ZggRealData,ka,Fac,'spline');
ImagZgg=interp2(Meshka,Meshga,ZggImagData,ka,Fac,'spline');

RealZbb=interp2(Meshka,Meshga,ZbbRealData,ka,Fac,'spline');
ImagZbb=interp2(Meshka,Meshga,ZbbImagData,ka,Fac,'spline');

RealZbt=interp2(Meshka,Meshga,ZbtRealData,ka,Fac,'spline');
ImagZbt=interp2(Meshka,Meshga,ZbtImagData,ka,Fac,'spline');

Zr1b = (RealZr1b + ImagZr1b*1i) * MultiFactor_L;
Zr1r1 = (RealZr1r1 + ImagZr1r1*1i) * MultiFactor_L;
Zr1g = (RealZr1g + ImagZr1g*1i) * MultiFactor_L;
Zr1r2 = (RealZr1r2 + ImagZr1r2*1i) * MultiFactor_L;
Zr1t = (RealZr1t + ImagZr1t*1i) * MultiFactor_L;

Zgb = (RealZgb + ImagZgb*1i) * MultiFactor_g;
Zgg = (RealZgg + ImagZgg*1i) * MultiFactor_g;

Zbb = (RealZbb + ImagZbb*1i) * MultiFactor_b;
Zbt = (RealZbt + ImagZbt*1i) * MultiFactor_b;

%%

% Construce Radiation Impednace Matrix
% {1,1} = r1_b, {1.2} = r1_r1, {1,3} = r1_g, {1,4} = r1_r2, {1,5} = r1_t
% {1,6} = g_b {1,7} = g_g1 {1,8} = b_b {1,9} = b_t

z_r=cell(1,length(ka)); % radiation impedance matrix

for i=1:length(ka);
    z_r{1,i}=zeros(zlength,zlength);
        
    %% 2Array
    % from 1st ring (Bottom Ring)
    z_r{1,i}(2,1)=-Zr1b(i)/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(2,2)=Zr1r1(i)/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(2,4)=Zr1g(i)/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(2,6)=Zr1r2(i)/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(2,7)=Zr1t(i)/(rhowater*sp*2*pi*a*L);
    
    %from 2nd ring (Top Ring)
    z_r{1,i}(6,1)=-z_r{1,i}(2,7);
    z_r{1,i}(6,2)=z_r{1,i}(2,6);
    z_r{1,i}(6,4)=z_r{1,i}(2,4);
    z_r{1,i}(6,6)=z_r{1,i}(2,2);
    z_r{1,i}(6,7)=-z_r{1,i}(2,1);
    
    %from gap
    z_r{1,i}(4,1)=-Zgb(i)/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(4,2)=Zr1g(i)/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(4,4)=Zgg(i)/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(4,6)=Zr1g(i)/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(4,7)=Zgb(i)/(rhowater*sp*2*pi*a*L);
    
    %from bottom, top
    z_r{1,i}(1,1)=Zbb(i)/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(1,2)=z_r{1,i}(2,1);
    z_r{1,i}(1,4)=z_r{1,i}(4,1);
    z_r{1,i}(1,6)=z_r{1,i}(6,1);
    z_r{1,i}(1,7)=-Zbt(i)/(rhowater*sp*2*pi*a*L);
    
    z_r{1,i}(7,1)=z_r{1,i}(1,7);
    z_r{1,i}(7,2)=z_r{1,i}(2,7);
    z_r{1,i}(7,4)=z_r{1,i}(4,7);
    z_r{1,i}(7,6)=z_r{1,i}(6,7);
    z_r{1,i}(7,7)=z_r{1,i}(1,1);
end
end