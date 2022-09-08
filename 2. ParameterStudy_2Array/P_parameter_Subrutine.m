function [pu,pp]=P_parameter_Subrutine(R0,Z0,Zam,ka,a,L,zlength2,AR,num)
% Define Normalized variables
eta0=R0/a; % Normalized radius
zta0=Z0/L; % Normalized FFR height
LaRatio = L/a;
% zta0=abs(Z0)/L; % Normalized height
% integral etas, ztas, theta

% Distance parameter
distance=@(etas,theta,ztas,LaRatio) sqrt(eta0^2+etas.^2-2.*eta0*etas.*cos(theta)+(LaRatio)^2.*(zta0-ztas).^2);
for Count = 1 : zlength2

if Count == 1
    z = Zam(1)/L;
    puintegral=@(etas,theta) (exp(-1i.*ka.*distance(etas,theta,z,LaRatio))./(distance(etas,theta,z,LaRatio))).*etas;
    pu(Count,1)=integral2(puintegral,0,1,0,2*pi,'Method','auto','AbsTol',0,'RelTol',1e-10).*((-1i*ka)/(4*pi));
    
    ppintegral=@(etas,theta) (((z-zta0).*(1+1i.*ka.*distance(etas,theta,z,LaRatio)))...
        .*((exp(-1i.*ka.*distance(etas,theta,z,LaRatio)))./((distance(etas,theta,z,LaRatio)).^3))).*etas;
    pp(Count,1)=integral2(ppintegral,0,1,0,2*pi,'Method','auto','AbsTol',0,'RelTol',1e-10).*(LaRatio/(4*pi));
elseif Count == zlength2
    z = Zam(end)/L;
    puintegral=@(etas,theta) (exp(-1i.*ka.*distance(etas,theta,z,LaRatio))./(distance(etas,theta,z,LaRatio))).*etas;
    pu(Count,1)=integral2(puintegral,0,1,0,2*pi,'Method','auto','AbsTol',0,'RelTol',1e-10).*((1i*ka)/(4*pi));
    
    ppintegral=@(etas,theta) (((z-zta0).*(1+1i.*ka.*distance(etas,theta,z,LaRatio)))...
        .*((exp(-1i.*ka.*distance(etas,theta,z,LaRatio)))./((distance(etas,theta,z,LaRatio)).^3))).*etas;
    pp(Count,1)=integral2(ppintegral,0,1,0,2*pi,'Method','auto','AbsTol',0,'RelTol',1e-10).*(-LaRatio/(4*pi));
else
    zmax = Zam(Count-1)/L;
    zmin = Zam(Count)/L;
%     zmax = Zam(Count)/L;
%     zmin = Zam(Count-1)/L;
    puintegral=@(ztas,theta) (exp(-1i.*ka.*distance(1,theta,ztas,LaRatio))./(distance(1,theta,ztas,LaRatio)));
    pu(Count,1)=integral2(puintegral,zmin,zmax,0,2*pi,'Method','auto','AbsTol',0,'RelTol',1e-10).*(((1i*ka)/(4*pi))*LaRatio);

    ppintegral=@(ztas,theta) (((1-eta0.*cos(theta)).*(1+1i.*ka.*distance(1,theta,ztas,LaRatio)))...
        .*((exp(-1i.*ka.*distance(1,theta,ztas,LaRatio)))./((distance(1,theta,ztas,LaRatio)).^3)));
    pp(Count,1)=integral2(ppintegral,zmin,zmax,0,2*pi,'Method','auto','AbsTol',0,'RelTol',1e-10).*(-LaRatio/(4*pi));
end

end

pp = (pp.'*AR).';

% Dimension expansion
switch num                  
    case 2        
       pu = [pu(1:2);0;pu(3);0;pu(4:5)];
       pp = [pp(1:2);0;pp(3);0;pp(4:5)];
    case 3
      pu = [pu(1:2);0;pu(3);0;pu(4);0;pu(5);0;pu(6:7)];
      pp = [pp(1:2);0;pp(3);0;pp(4);0;pp(5);0;pp(6:7)];
    case 4
      pu = [pu(1:2);0;pu(3);0;pu(4);0;pu(5);0;pu(6);0;pu(7);0;pu(8:9)];
      pp = [pp(1:2);0;pp(3);0;pp(4);0;pp(5);0;pp(6);0;pp(7);0;pp(8:9)];
end   
