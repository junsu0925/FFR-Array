function [Vu_Dimensionless]=HKI_Sub_VelVec(VelSurSel, u0, NN)
Nb = NN(1); Nr = NN(2); Nt = NN(3);

if VelSurSel == 1;
    % Bottom 면만 1m/s로 구동
    Vu = zeros(Nb+Nr+Nt,1);
    Vu(1:Nb,1) = -1;
elseif VelSurSel == 2;
    % Cylinder 면만 1m/s로 구동
    Vu = zeros(Nb+Nr+Nt,1);
    Vu(Nb+1:Nb+Nr,1) = 1;
elseif VelSurSel == 3;
    % Top 면만 1m/s로 구동
    Vu = zeros(Nb+Nr+Nt,1);
    Vu(Nb+Nr+1:end,1) = 1;
elseif VelSurSel == 4;
    % 모든 surface 가 1m/s 로 구동
    Vu = ones(Nb+Nr+Nt,1);
    Vu(1:Nb,1) = -1;
end

%% Dimensionless
Vu_Dimensionless = Vu./u0;
end