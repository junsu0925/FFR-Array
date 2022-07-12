function [z_r,exceldata]=Radiation_Impedance_3Array_Subrutine(rhowater,sp,ka,a,L,zlength)
% Calculate z_sR [Page 48]

% Import Data from Excel File
% sheet numbers -> %

%% 3-Array

exceldata=cell(1,19);
for i=1:19;
    [~, exceldata{1,i}]=xlsread('3array.xlsx',i,'B6:B801');
    for j=1:length(ka);
        exceldata{1,i}{j,1}=str2double(exceldata{1,i}{j,1});
    end
end

%% Construce Radiation Impednace Matrix
z_r=cell(1,length(ka)); % radiation impedance matrix

for i=1:length(ka);
    z_r{1,i}=zeros(zlength,zlength);

    %% 3Array
    % from Top
    z_r{1,i}(1,1)=exceldata{1,1}{i,1}/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(1,2)=exceldata{1,2}{i,1}/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(1,4)=exceldata{1,3}{i,1}/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(1,6)=exceldata{1,4}{i,1}/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(1,8)=exceldata{1,5}{i,1}/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(1,10)=exceldata{1,6}{i,1}/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(1,11)=-exceldata{1,7}{i,1}/(rhowater*sp*2*pi*a*L);
    
    % from 1st ring (Top ring)
    z_r{1,i}(2,1)=z_r{1,i}(1,2);
    z_r{1,i}(2,2)=exceldata{1,8}{i,1}/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(2,4)=exceldata{1,9}{i,1}/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(2,6)=exceldata{1,10}{i,1}/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(2,8)=exceldata{1,11}{i,1}/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(2,10)=exceldata{1,12}{i,1}/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(2,11)=-exceldata{1,13}{i,1}/(rhowater*sp*2*pi*a*L);
    
    % from 1st gab
    z_r{1,i}(4,1)=z_r{1,i}(1,4);
    z_r{1,i}(4,2)=z_r{1,i}(2,4);
    z_r{1,i}(4,4)=exceldata{1,14}{i,1}/(rhowater*sp*2*pi*a*L); 
    z_r{1,i}(4,6)=exceldata{1,15}{i,1}/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(4,8)=exceldata{1,16}{i,1}/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(4,10)=z_r{1,i}(2,8);
    z_r{1,i}(4,11)=-z_r{1,i}(1,8);
    
    % from 2nd ring (Mid ring)
    z_r{1,i}(6,1)=z_r{1,i}(1,6);
    z_r{1,i}(6,2)=z_r{1,i}(2,6);
    z_r{1,i}(6,4)=z_r{1,i}(4,6); 
    z_r{1,i}(6,6)=exceldata{1,17}{i,1}/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(6,8)=z_r{1,i}(4,6);
    z_r{1,i}(6,10)=z_r{1,i}(2,6);
    z_r{1,i}(6,11)=-z_r{1,i}(1,6);
    
    % from 2nd gab
    z_r{1,i}(8,1)=z_r{1,i}(1,8);
    z_r{1,i}(8,2)=z_r{1,i}(2,8);
    z_r{1,i}(8,4)=z_r{1,i}(4,8); 
    z_r{1,i}(8,6)=z_r{1,i}(4,6);
    z_r{1,i}(8,8)=exceldata{1,18}{i,1}/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(8,10)=z_r{1,i}(2,4);
    z_r{1,i}(8,11)=-z_r{1,i}(1,4);
    
    % from 3nd ring (Bottom ring)
    z_r{1,i}(10,1)=z_r{1,i}(1,10);
    z_r{1,i}(10,2)=z_r{1,i}(2,10);
    z_r{1,i}(10,4)=z_r{1,i}(2,8); 
    z_r{1,i}(10,6)=z_r{1,i}(2,6);
    z_r{1,i}(10,8)=z_r{1,i}(2,4);
    z_r{1,i}(10,10)=exceldata{1,19}{i,1}/(rhowater*sp*2*pi*a*L);
    z_r{1,i}(10,11)=-z_r{1,i}(1,2);

    % from bottom 
    z_r{1,i}(11,1)=z_r{1,i}(1,11);
    z_r{1,i}(11,2)=z_r{1,i}(2,11);
    z_r{1,i}(11,4)=z_r{1,i}(4,11); 
    z_r{1,i}(11,6)=z_r{1,i}(6,11);
    z_r{1,i}(11,8)=z_r{1,i}(8,11);
    z_r{1,i}(11,10)=z_r{1,i}(10,11);
    z_r{1,i}(11,11)=z_r{1,i}(1,1);
        
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