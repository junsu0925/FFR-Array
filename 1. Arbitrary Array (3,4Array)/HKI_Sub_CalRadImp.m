function [VF_For_FFR] = HKI_Sub_CalRadImp(InvMat_HKI,Geometry,NN,Num)

radius_a = Geometry.a;
length_l = Geometry.RL;
gap_g = Geometry.g(1);

Nb = NN(1); Nt = NN(end);

TotNumPre = sum(NN)-2*Num;
TotNumVel = sum(NN);

%% Radiation Impedance Matrix Define (same as Force)

T_3B = zeros(TotNumVel,2*Num+1);

for i = 1 : 2*Num+1
    if i == 1
        Sum_Num_start = 1;
        Sum_Num_end = NN(1);
    else
        Sum_Num_start = 1+sum(NN(1:i-1));
        Sum_Num_end = sum(NN(1:i));
    end
    T_3B(Sum_Num_start:Sum_Num_end,i) = ones(NN(i),1);        
end

GA = zeros(2*Num+1,TotNumPre);
VA_b = zeros(1,Nb);
% VA_r = zeros(1,NN(2));
VA_t = zeros(1,NN(end));
%% Area Vector define
VA_b(1,1) = 1/6;
VA_b(1,2:Nb-1) = (1:Nb-2);
VA_b(1,end) = (3*(Nb-1)-1)/6;
VA_b = -VA_b.*((2*pi)/(Nb-1)^2);

VA_t(1,1) =(3*(Nt-1)-1)/6;
VA_t(1,2:Nt-1) = (1:Nt-2);
VA_t(1,end) = 1/6;
VA_t = VA_t.*((2*pi)/(Nt-1)^2);


for i = 1 : 2*Num+1
    if i == 1
        GA(i,1:NN(1)) = VA_b;
    elseif i == 2*Num+1
        Sum_Num_start = sum(NN(1:i-1))-i+2;
        Sum_Num_end = sum(NN(1:i))-i+1;

        GA(i,Sum_Num_start:Sum_Num_end) = VA_t;
    else
        Sum_Num_start = sum(NN(1:i-1))-i+2;
        Sum_Num_end = sum(NN(1:i))-i+1;

        VA_r = zeros(1,NN(i));
        VA_r(1,1) = 1/2;
        VA_r(1,2:NN(i)-1) = 1;
        VA_r(1,end) = 1/2;
        VA_r = VA_r.*((2*pi)/(NN(i)-1));

        GA(i,Sum_Num_start:Sum_Num_end) = VA_r;      
    end    
end

%% Radiation Impedance Matrix Define (same as Force)
HKI_to_FFR_AreaRatio_M = zeros(2*Num+1,2*Num+1);
switch Num
    case 1
        HKI_to_FFR_AreaRatio_M(1,1) = radius_a/(2*pi*length_l);
        HKI_to_FFR_AreaRatio_M(2,2) = 1/(2*pi);
        HKI_to_FFR_AreaRatio_M(3,3) = radius_a/(2*pi*length_l);                
    case 2        
        HKI_to_FFR_AreaRatio_M(1,1) = radius_a/(2*pi*length_l);
        HKI_to_FFR_AreaRatio_M(2,2) = 1/(2*pi);
        HKI_to_FFR_AreaRatio_M(3,3) = gap_g/(2*pi*length_l);
        HKI_to_FFR_AreaRatio_M(4,4) = 1/(2*pi);
        HKI_to_FFR_AreaRatio_M(5,5) = radius_a/(2*pi*length_l);  
    case 3
        HKI_to_FFR_AreaRatio_M(1,1) = radius_a/(2*pi*length_l);
        HKI_to_FFR_AreaRatio_M(2,2) = 1/(2*pi);
        HKI_to_FFR_AreaRatio_M(3,3) = gap_g/(2*pi*length_l);
        HKI_to_FFR_AreaRatio_M(4,4) = 1/(2*pi);
        HKI_to_FFR_AreaRatio_M(5,5) = gap_g/(2*pi*length_l);
        HKI_to_FFR_AreaRatio_M(6,6) = 1/(2*pi);
        HKI_to_FFR_AreaRatio_M(7,7) = radius_a/(2*pi*length_l);
    case 4
        HKI_to_FFR_AreaRatio_M(1,1) = radius_a/(2*pi*length_l);
        HKI_to_FFR_AreaRatio_M(2,2) = 1/(2*pi);
        HKI_to_FFR_AreaRatio_M(3,3) = gap_g/(2*pi*length_l);
        HKI_to_FFR_AreaRatio_M(4,4) = 1/(2*pi);
        HKI_to_FFR_AreaRatio_M(5,5) = gap_g/(2*pi*length_l);
        HKI_to_FFR_AreaRatio_M(6,6) = 1/(2*pi);
        HKI_to_FFR_AreaRatio_M(7,7) = gap_g/(2*pi*length_l);
        HKI_to_FFR_AreaRatio_M(8,8) = 1/(2*pi);
        HKI_to_FFR_AreaRatio_M(9,9) = radius_a/(2*pi*length_l); 
end

VF_For_FFR = HKI_to_FFR_AreaRatio_M*GA*InvMat_HKI*T_3B;

%% Area notation issue (HKI & FFR)
VF_For_FFR(:,1) = -VF_For_FFR(:,1);
end


