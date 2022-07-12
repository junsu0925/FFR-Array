function [z_r,exceldata]=Radiation_Impedance_Single_Subrutine(rhowater,sp,ka,a,L,zlength)
% Calculate z_sR [Page 48]

% Import Data from Excel File
% sheet numbers -> %
exceldata=cell(1,5);
for i=1:5;
    [~, exceldata{1,i}]=xlsread('Single.xlsx',i,'B6:B3981');
    for j=1:length(ka);
        exceldata{1,i}{j,1}=str2double(exceldata{1,i}{j,1});
    end
end

% Construce Radiation Impednace Matrix
z_r=cell(1,length(ka)); % radiation impedance matrix

for i=1:length(ka);
    z_r{1,i}=zeros(zlength,zlength);
   
    %% Original 
    % from 1st ring
    z_r{1,i}(2,1)=exceldata{1,1}{i,1}/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(2,2)=exceldata{1,2}{i,1}/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(2,3)=-exceldata{1,1}{i,1}/(rhowater*sp*2*pi*a*L);

    %from bottom, top
    z_r{1,i}(1,1)=exceldata{1,3}{i,1}/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(1,2)=exceldata{1,1}{i,1}/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(1,3)=-exceldata{1,5}{i,1}/(rhowater*sp*2*pi*a*L);
    
    z_r{1,i}(3,1)=z_r{1,i}(1,3);
    z_r{1,i}(3,2)=-exceldata{1,1}{i,1}/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(3,3)=z_r{1,i}(1,1);   
    
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