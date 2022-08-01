function [z_rr,z_rz,z_1,z_2,z_PRm]=FFR_Sub_CavImp(LaRatio,t,a,ka,rhopzt,rhowater,sp,sE11,n,rootj0dot)
% 변수 변환
% f -> ka
% a, L -> LaRatio
% t, a -> Ratio

% omega=2*pi*f;
% k=@(w) w/sp;

% Calculate z_rr
z_rrpart=@(nka, nn) ka(nka).^2./((ka(nka).^2-(rootj0dot(nn)).^2));
z_rr=zeros(1,length(ka));
z_rrtemp=zeros(1,length(ka));

for i=1:n;
    for j=1:length(ka);
        z_rrtemp(j)=z_rrpart(j,i);
        z_rr(j)=z_rr(j)+z_rrtemp(j);
    end
end

clear z_rrtemp;

for i=1:length(ka);
    z_rr(i)=(-1i*2/ka(i))*(1+z_rr(i));
end

% Caculate z_rz
z_rz=zeros(1,length(ka));
for i=1:length(ka);
    z_rz(i)=(1/(1i*ka(i)))*(1/LaRatio);
end

% Caculate z_1
z_1part=@(nka,nn) ka(nka).^2./(ka(nka).^2-(((nn*pi).^2)./LaRatio^2));
z_1=zeros(1,length(ka));
z_1temp=zeros(1,length(ka));

for i=1:n;
    for j=1:length(ka);
        z_1temp(j)=z_1part(j,i);
        z_1(j)=z_1(j)+z_1temp(j);
    end
end

clear Z1temp;

for i=1:length(ka);
    z_1(i)=(-1i/(2*ka(i)*LaRatio^2))*(1+2*z_1(i));
end

% Caculate z_2
% z_2part=@(nka,nn) (((-1).^nn).*ka(nka).^2)./(ka(nka).^2-((nn*pi).^2./LaRatio)^2);
% 20161202 수정
z_2part=@(nka,nn) (((-1).^nn).*ka(nka).^2)./(ka(nka).^2-((nn*pi)./LaRatio)^2);
z_2=zeros(1,length(ka));
z_2temp=zeros(1,length(ka));

for i=1:n;
    for j=1:length(ka);
        z_2temp(j)=z_2part(j,i);
        z_2(j)=z_2(j)+z_2temp(j);
    end
end

clear Z2temp;

for i=1:length(ka);
    z_2(i)=(-1i./(2.*ka(i)*LaRatio^2)).*(1+2.*z_2(i));
end

% Calculate z_PR
Ratio = t/a;
z_PRmpart=@(nka) (1i.*ka(nka).*(Ratio)*(rhopzt/rhowater))+(Ratio/(1i.*ka(nka).*(sE11*rhowater*sp^2)));
z_PRm=zeros(1,length(ka));

for i=1:length(ka);
    z_PRm(i)=z_PRmpart(i);
end
end