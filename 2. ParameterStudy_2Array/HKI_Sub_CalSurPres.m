function [Mat_3] = HKI_Sub_CalSurPres(GD, GW)

%% Define Velocity Vector
IdentityGW = eye(length(GW));
Mat_1 = (IdentityGW - GW);
Mat_3 = Mat_1\GD;

end