function [p_1u, p_r1u, p_r2u, p_2u, p_gu, p_1p, p_r1p, p_r2p, p_gp, p_2p]=Pres_para_2Array(R0,Z0,ka,LaRatio,a,L,g)
% 2Array Case
% Define Normalized variables
eta0=R0/a; % Normalized radius
zta0=Z0/L; % Normalized FFR height
gLRatio=g/L; % Normalized Gab and FFR hight
gaRatio=g/a; % Normalized Gab higjt and Gab radius
% zta0=abs(Z0)/L; % Normalized height
% integral etas, ztas, theta

% Distance parameter
distance=@(etas,theta,ztas,LaRatio) sqrt(eta0^2+etas.^2-2.*eta0*etas.*cos(theta)+(LaRatio)^2.*(zta0-ztas).^2);

%%
% Calculate p_1u
p_1u=zeros(1,length(ka));
zmin = -(1/2)* gLRatio - 1;
for i=1:length(ka);
%     p_1uintegral=@(etas,theta) (exp(-1i*ka(i).*distance_p1(etas,theta))./distance_p1(etas,theta)).*etas;
    p_1uintegral=@(etas,theta) (exp(-1i.*ka(i).*distance(etas,theta,zmin,LaRatio))./(distance(etas,theta,zmin,LaRatio))).*etas;
    p_1u(i)=integral2(p_1uintegral,0,1,0,2*pi,'Method','auto','AbsTol',0,'RelTol',1e-10).*((-1i*ka(i))/(4*pi));
end

% Calculate p_r1u
p_r1u=zeros(1,length(ka));
zmin = -(1/2)* gLRatio - 1;
zmax = -(1/2)* gLRatio;
for i=1:length(ka);
    p_r1uintegral=@(ztas,theta) (exp(-1i.*ka(i).*distance(1,theta,ztas,LaRatio))./(distance(1,theta,ztas,LaRatio)));
    p_r1u(i)=integral2(p_r1uintegral,zmin,zmax,0,2*pi,'Method','auto','AbsTol',0,'RelTol',1e-10).*(((1i*ka(i))/(4*pi))*LaRatio);
end

% Calculate p_r2u
p_r2u=zeros(1,length(ka));
zmin = (1/2)* gLRatio;
zmax = (1/2)* gLRatio + 1;
for i=1:length(ka);
    p_r2uintegral=@(ztas,theta) (exp(-1i.*ka(i).*distance(1,theta,ztas,LaRatio))./(distance(1,theta,ztas,LaRatio)));
    p_r2u(i)=integral2(p_r2uintegral,zmin,zmax,0,2*pi,'Method','auto','AbsTol',0,'RelTol',1e-10).*(((1i*ka(i))/(4*pi))*LaRatio);
end

% Calculate p_gu
p_gu=zeros(1,length(ka));
zmin = -(1/2)* gLRatio;
zmax = (1/2)* gLRatio;
for i=1:length(ka);
    p_guintegral=@(ztas,theta) (exp(-1i.*ka(i).*distance(1,theta,ztas,LaRatio))./(distance(1,theta,ztas,LaRatio)));
    p_gu(i)=integral2(p_guintegral,zmin,zmax,0,2*pi,'Method','auto','AbsTol',0,'RelTol',1e-10).*(((1i*ka(i))/(4*pi))*LaRatio);
end

% for i=1:length(ka);
%     p_guintegral=@(ztas,theta) (exp(-1i.*ka(i).*distance(1,theta,ztas,gaRatio))./(distance(1,theta,ztas,gaRatio)));
%     p_gu(i)=integral2(p_guintegral,zmin,zmax,0,2*pi,'Method','auto','AbsTol',0,'RelTol',1e-10).*(((1i*ka(i))/(4*pi))*gaRatio);
% end

% Calculate p_2u
p_2u=zeros(1,length(ka));
zmax = (1/2)* gLRatio + 1;
for i=1:length(ka);
    p_2uintegral=@(etas,theta) (exp(-1i.*ka(i).*distance(etas,theta,zmax,LaRatio))./(distance(etas,theta,zmax,LaRatio))).*etas;
    p_2u(i)=integral2(p_2uintegral,0,1,0,2*pi,'Method','auto','AbsTol',0,'RelTol',1e-10).*((1i*ka(i))/(4*pi));
end

%%
% Calculate p_1p
p_1p=zeros(1,length(ka));
zmin = -(1/2)* gLRatio - 1;
for i=1:length(ka);
    p_1pintegral=@(etas,theta) (((zmin-zta0).*(1+1i.*ka(i).*distance(etas,theta,zmin,LaRatio)))...
        .*((exp(-1i.*ka(i).*distance(etas,theta,-1/2,LaRatio)))./((distance(etas,theta,zmin,LaRatio)).^3))).*etas;
%     p_1pintegral=@(etas,theta) ((((1/2)-eta0).*(1+1i.*ka(i).*distance(etas,theta,-1/2,LaRatio)))...
%         .*((exp(-1i.*ka(i).*distance(etas,theta,-1/2,LaRatio)))./((distance(etas,theta,-1/2,LaRatio)).^3))).*etas;
    p_1p(i)=integral2(p_1pintegral,0,1,0,2*pi,'Method','auto','AbsTol',0,'RelTol',1e-10).*(LaRatio/(4*pi));
end

% Calculate p_r1p
p_r1p=zeros(1,length(ka));
zmin = -(1/2)* gLRatio - 1;
zmax = -(1/2)* gLRatio;
for i=1:length(ka);
    p_r1pintegral=@(ztas,theta) (((1-eta0.*cos(theta)).*(1+1i.*ka(i).*distance(1,theta,ztas,LaRatio)))...
        .*((exp(-1i.*ka(i).*distance(1,theta,ztas,LaRatio)))./((distance(1,theta,ztas,LaRatio)).^3)));
    p_r1p(i)=integral2(p_r1pintegral,zmin,zmax,0,2*pi,'Method','auto','AbsTol',0,'RelTol',1e-10).*(-LaRatio/(4*pi));
end

% Calculate p_r2p
p_r2p=zeros(1,length(ka));
zmin = (1/2)* gLRatio;
zmax = (1/2)* gLRatio + 1;
for i=1:length(ka);
    p_r2pintegral=@(ztas,theta) (((1-eta0.*cos(theta)).*(1+1i.*ka(i).*distance(1,theta,ztas,LaRatio)))...
        .*((exp(-1i.*ka(i).*distance(1,theta,ztas,LaRatio)))./((distance(1,theta,ztas,LaRatio)).^3)));
    p_r2p(i)=integral2(p_r2pintegral,zmin,zmax,0,2*pi,'Method','auto','AbsTol',0,'RelTol',1e-10).*(-LaRatio/(4*pi));
end

% Calculate p_gp
p_gp=zeros(1,length(ka));
zmin = -(1/2)* gLRatio;
zmax = (1/2)* gLRatio;
for i=1:length(ka);
    p_gpintegral=@(ztas,theta) (((1-eta0.*cos(theta)).*(1+1i.*ka(i).*distance(1,theta,ztas,LaRatio)))...
        .*((exp(-1i.*ka(i).*distance(1,theta,ztas,LaRatio)))./((distance(1,theta,ztas,LaRatio)).^3)));
    p_gp(i)=integral2(p_gpintegral,zmin,zmax,0,2*pi,'Method','auto','AbsTol',0,'RelTol',1e-10).*(-LaRatio/(4*pi));
end

% for i=1:length(ka);
%     p_gpintegral=@(ztas,theta) (((1-eta0.*cos(theta)).*(1+1i.*ka(i).*distance(1,theta,ztas,LaRatio)))...
%         .*((exp(-1i.*ka(i).*distance(1,theta,ztas,LaRatio)))./((distance(1,theta,ztas,LaRatio)).^3)));
%     p_gp(i)=integral2(p_gpintegral,zmin,zmax,0,2*pi,'Method','auto','AbsTol',0,'RelTol',1e-10).*(-gaRatio/(4*pi));
% end

% for i=1:length(ka);
%     p_gpintegral=@(ztas,theta) (((1-eta0.*cos(theta)).*(1+1i.*ka(i).*distance(1,theta,ztas,gaRatio)))...
%         .*((exp(-1i.*ka(i).*distance(1,theta,ztas,gaRatio)))./((distance(1,theta,ztas,gaRatio)).^3)));
%     p_gp(i)=integral2(p_gpintegral,zmin,zmax,0,2*pi,'Method','auto','AbsTol',0,'RelTol',1e-10).*(-gaRatio/(4*pi));
% end

% Calculate p_2p
p_2p=zeros(1,length(ka));
zmax = (1/2)* gLRatio + 1;
for i=1:length(ka);
    p_2pintegral=@(etas,theta) (((zmax-zta0).*(1+1i.*ka(i).*distance(etas,theta,zmax,LaRatio)))...
        .*((exp(-1i.*ka(i).*distance(etas,theta,1/2,LaRatio)))./((distance(etas,theta,zmax,LaRatio)).^3))).*etas;
    p_2p(i)=integral2(p_2pintegral,0,1,0,2*pi,'Method','auto','AbsTol',0,'RelTol',1e-10).*(-LaRatio/(4*pi));
end
end