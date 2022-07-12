function [z_r,exceldata]=Radiation_Impedance_4Array_Subrutine(rhowater,sp,ka,a,L,zlength)
% Calculate z_sR [Page 48]

% Import Data from Excel File
% sheet numbers -> %

%% 4-Array

exceldata=cell(1,28);
for i=1:28;
    [~, exceldata{1,i}]=xlsread('4array.xlsx',i,'B6:B801');
    for j=1:length(ka);
        exceldata{1,i}{j,1}=str2double(exceldata{1,i}{j,1});
    end
end

%% Construce Radiation Impednace Matrix
z_r=cell(1,length(ka)); % radiation impedance matrix

for i=1:length(ka);
    z_r{1,i}=zeros(zlength,zlength);

    %% 4Array
    % from Top
    z_r{1,i}(1,1)=exceldata{1,1}{i,1}/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(1,2)=exceldata{1,2}{i,1}/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(1,4)=exceldata{1,3}{i,1}/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(1,6)=exceldata{1,4}{i,1}/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(1,8)=exceldata{1,5}{i,1}/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(1,10)=exceldata{1,6}{i,1}/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(1,12)=exceldata{1,7}{i,1}/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(1,14)=exceldata{1,8}{i,1}/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(1,15)=-exceldata{1,9}{i,1}/(rhowater*sp*2*pi*a*L);
    
    % from 1st ring (Top ring)
    z_r{1,i}(2,1)=z_r{1,i}(1,2);
    z_r{1,i}(2,2)=exceldata{1,10}{i,1}/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(2,4)=exceldata{1,11}{i,1}/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(2,6)=exceldata{1,12}{i,1}/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(2,8)=exceldata{1,13}{i,1}/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(2,10)=exceldata{1,14}{i,1}/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(2,12)=exceldata{1,15}{i,1}/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(2,14)=exceldata{1,16}{i,1}/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(2,15)=-z_r{1,i}(1,14);
    
    % from 1st gab
    z_r{1,i}(4,1)=z_r{1,i}(1,4);
    z_r{1,i}(4,2)=z_r{1,i}(2,4);
    z_r{1,i}(4,4)=exceldata{1,17}{i,1}/(rhowater*sp*2*pi*a*L); 
    z_r{1,i}(4,6)=exceldata{1,18}{i,1}/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(4,8)=exceldata{1,19}{i,1}/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(4,10)=exceldata{1,20}{i,1}/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(4,12)=exceldata{1,21}{i,1}/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(4,14)=z_r{1,i}(2,12);
    z_r{1,i}(4,15)=-z_r{1,i}(1,12);
    
    % from 2nd ring (Mid 1 ring)
    z_r{1,i}(6,1)=z_r{1,i}(1,6);
    z_r{1,i}(6,2)=z_r{1,i}(2,6);
    z_r{1,i}(6,4)=z_r{1,i}(4,6); 
    z_r{1,i}(6,6)=exceldata{1,22}{i,1}/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(6,8)=exceldata{1,23}{i,1}/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(6,10)=exceldata{1,24}{i,1}/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(6,12)=z_r{1,i}(4,10);
    z_r{1,i}(6,14)=z_r{1,i}(2,10);
    z_r{1,i}(6,15)=-z_r{1,i}(1,10);
    
    % from 2nd gab
    z_r{1,i}(8,1)=z_r{1,i}(1,8);
    z_r{1,i}(8,2)=z_r{1,i}(2,8);
    z_r{1,i}(8,4)=z_r{1,i}(4,8); 
    z_r{1,i}(8,6)=z_r{1,i}(6,8);
    z_r{1,i}(8,8)=exceldata{1,25}{i,1}/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(8,10)=z_r{1,i}(6,8);
    z_r{1,i}(8,12)=z_r{1,i}(4,8);
    z_r{1,i}(8,14)=z_r{1,i}(2,8);
    z_r{1,i}(8,15)=-z_r{1,i}(1,8);
    
    % from 3nd ring (Mid 2 ring)
    z_r{1,i}(10,1)=z_r{1,i}(1,10);
    z_r{1,i}(10,2)=z_r{1,i}(2,10);
    z_r{1,i}(10,4)=z_r{1,i}(4,10); 
    z_r{1,i}(10,6)=z_r{1,i}(6,10);
    z_r{1,i}(10,8)=z_r{1,i}(8,10);
    z_r{1,i}(10,10)=exceldata{1,26}{i,1}/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(10,12)=z_r{1,i}(4,6);
    z_r{1,i}(10,14)=z_r{1,i}(2,6);
    z_r{1,i}(10,15)=-z_r{1,i}(1,6);

    % from 3rd gab 
    z_r{1,i}(12,1)=z_r{1,i}(1,12);
    z_r{1,i}(12,2)=z_r{1,i}(2,12);
    z_r{1,i}(12,4)=z_r{1,i}(4,12); 
    z_r{1,i}(12,6)=z_r{1,i}(6,12);
    z_r{1,i}(12,8)=z_r{1,i}(8,12);
    z_r{1,i}(12,10)=z_r{1,i}(10,12);
    z_r{1,i}(12,12)=exceldata{1,27}{i,1}/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(12,14)=z_r{1,i}(2,4);
    z_r{1,i}(12,15)=-z_r{1,i}(1,4);
    
    % from 4th ring
    z_r{1,i}(14,1)=z_r{1,i}(1,14);
    z_r{1,i}(14,2)=z_r{1,i}(2,14);
    z_r{1,i}(14,4)=z_r{1,i}(4,14); 
    z_r{1,i}(14,6)=z_r{1,i}(6,14);
    z_r{1,i}(14,8)=z_r{1,i}(8,14);
    z_r{1,i}(14,10)=z_r{1,i}(10,14);
    z_r{1,i}(14,12)=z_r{1,i}(12,14);
    z_r{1,i}(14,14)=exceldata{1,28}{i,1}/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(14,15)=-z_r{1,i}(1,2);
    
    % from bottom
    z_r{1,i}(15,1)=z_r{1,i}(1,15);
    z_r{1,i}(15,2)=z_r{1,i}(2,15);
    z_r{1,i}(15,4)=z_r{1,i}(4,15); 
    z_r{1,i}(15,6)=z_r{1,i}(6,15);
    z_r{1,i}(15,8)=z_r{1,i}(8,15);
    z_r{1,i}(15,10)=z_r{1,i}(10,15);
    z_r{1,i}(15,12)=z_r{1,i}(12,15);
    z_r{1,i}(15,14)=z_r{1,i}(14,15);
    z_r{1,i}(15,15)=z_r{1,i}(1,1);
        
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