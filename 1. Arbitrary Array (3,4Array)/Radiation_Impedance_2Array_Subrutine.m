function [z_r,exceldata]=Radiation_Impedance_2Array_Subrutine(rhowater,sp,ka,a,L,zlength)
% Calculate z_sR [Page 48]

% Import Data from Excel File
% sheet numbers -> %

%% 2-Array

exceldata=cell(1,9);
for i=1:9;
    [~, exceldata{1,i}]=xlsread('2array.xlsx',i,'B6:B801');
    for j=1:length(ka);
        exceldata{1,i}{j,1}=str2double(exceldata{1,i}{j,1});
    end
end

%% Construce Radiation Impednace Matrix
z_r=cell(1,length(ka)); % radiation impedance matrix

for i=1:length(ka);
    z_r{1,i}=zeros(zlength,zlength);
       
    %% 2Array  
    % from 1st ring (Top Ring)
    z_r{1,i}(2,1)=exceldata{1,1}{i,1}/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(2,2)=exceldata{1,2}{i,1}/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(2,4)=exceldata{1,3}{i,1}/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(2,6)=exceldata{1,4}{i,1}/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(2,7)=-exceldata{1,5}{i,1}/(rhowater*sp*2*pi*a*L);
    
    %from 2nd ring (Bottom Ring)
    z_r{1,i}(6,1)=-z_r{1,i}(2,7);
    z_r{1,i}(6,2)=z_r{1,i}(2,6);
    z_r{1,i}(6,4)=z_r{1,i}(2,4);
    z_r{1,i}(6,6)=z_r{1,i}(2,2);
    z_r{1,i}(6,7)=-z_r{1,i}(2,1);
    
    %from gap
    z_r{1,i}(4,1)=exceldata{1,6}{i,1}/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(4,2)=exceldata{1,3}{i,1}/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(4,4)=exceldata{1,7}{i,1}/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(4,6)=exceldata{1,3}{i,1}/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(4,7)=-exceldata{1,6}{i,1}/(rhowater*sp*2*pi*a*L);
    
    %from Top, Bottom
    z_r{1,i}(1,1)=exceldata{1,8}{i,1}/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(1,2)=z_r{1,i}(2,1);
    z_r{1,i}(1,4)=z_r{1,i}(4,1);
    z_r{1,i}(1,6)=z_r{1,i}(6,1);
    z_r{1,i}(1,7)=-exceldata{1,9}{i,1}/(rhowater*sp*2*pi*a*L);
    
    z_r{1,i}(7,1)=z_r{1,i}(1,7);
    z_r{1,i}(7,2)=z_r{1,i}(2,7);
    z_r{1,i}(7,4)=z_r{1,i}(4,7);
    z_r{1,i}(7,6)=z_r{1,i}(6,7);
    z_r{1,i}(7,7)=z_r{1,i}(1,1);
        
        
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