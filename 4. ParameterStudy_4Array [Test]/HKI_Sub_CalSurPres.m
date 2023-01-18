function [Mat_3] = HKI_Sub_CalSurPres(GD, GW,etaOb,ztaOb)

%% Define Velocity Vector
IdentityGW = eye(length(GW)-length(etaOb));
Modified_IdentityGW = [IdentityGW ; zeros(length(etaOb),length(IdentityGW))];
Mat_1 = (Modified_IdentityGW - GW);
Mat_3 = Mat_1\GD;

end