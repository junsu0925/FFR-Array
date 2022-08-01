function [z_r,VF_For_FFR]=FFR_Sub_RadImp(VF_For_FFR,rhowater,sp,ka,a,L,zlength)
% Calculate z_sR [Page 48]

% Import Data from Excel File
% sheet numbers -> %
% exceldata=cell(1,5);
% for i=1:5;
%     [~, exceldata{1,i}]=xlsread('Single.xlsx',i,'B6:B3006');
%     for j=1:length(ka);
%         exceldata{1,i}{j,1}=str2double(exceldata{1,i}{j,1});
%     end
% end

HKI_data=cell(1,5);
HKI_data{1,1} = VF_For_FFR(:,1);
HKI_data{1,2} = VF_For_FFR(:,2);
HKI_data{1,3} = VF_For_FFR(:,3);
HKI_data{1,4} = VF_For_FFR(:,4);
HKI_data{1,5} = VF_For_FFR(:,5);

% Construce Radiation Impednace Matrix
z_r=cell(1,length(ka)); % radiation impedance matrix

for i=1:length(ka);
    z_r{1,i}=zeros(zlength,zlength);
    
    %% Original 
    % from 1st ring
    z_r{1,i}(2,1)=-HKI_data{1,1}(i,1);
    z_r{1,i}(2,2)=HKI_data{1,2}(i,1);
    z_r{1,i}(2,3)=HKI_data{1,1}(i,1);

    %from bottom, top
    z_r{1,i}(1,1)=HKI_data{1,3}(i,1);
    z_r{1,i}(1,2)=-HKI_data{1,1}(i,1);
    z_r{1,i}(1,3)=-HKI_data{1,5}(i,1);
    
    z_r{1,i}(3,1)=z_r{1,i}(1,3);
    z_r{1,i}(3,2)=HKI_data{1,1}(i,1);
    z_r{1,i}(3,3)=z_r{1,i}(1,1); 
   
%     %% Original 
%     % from 1st ring
%     z_r{1,i}(2,1)=-exceldata{1,1}{i,1}/(rhowater*sp*2*pi*a*L);
%     z_r{1,i}(2,2)=exceldata{1,2}{i,1}/(rhowater*sp*2*pi*a*L);
%     z_r{1,i}(2,3)=exceldata{1,1}{i,1}/(rhowater*sp*2*pi*a*L);
% 
%     %from bottom, top
%     z_r{1,i}(1,1)=exceldata{1,3}{i,1}/(rhowater*sp*2*pi*a*L);
%     z_r{1,i}(1,2)=-exceldata{1,1}{i,1}/(rhowater*sp*2*pi*a*L);
%     z_r{1,i}(1,3)=-exceldata{1,5}{i,1}/(rhowater*sp*2*pi*a*L);
%     
%     z_r{1,i}(3,1)=z_r{1,i}(1,3);
%     z_r{1,i}(3,2)=exceldata{1,1}{i,1}/(rhowater*sp*2*pi*a*L);
%     z_r{1,i}(3,3)=z_r{1,i}(1,1);   
    
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