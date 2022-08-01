function [R]=HKI_Sub_R(eta_o,zta_o,eta,theta,zta,LaRatio)
R = sqrt(eta_o^2 + eta.^2 - 2*eta_o.*eta.*cos(theta) + (LaRatio)^2.*(zta_o-zta).^2);
end