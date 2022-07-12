function [RadiationTable,RadiationAreaFactorTable,TableParameter] = RadiationImpTable_ParameterNo5(Geometry, rhowater, sp, kaForTable)
% {1,1} = t_t, {1.2} = t_r1, {1,3} = t_g1, {1,4} = t_r2, {1,5} = t_g2,
% {1,6} = t_r3, {1,7} = t_g3, {1,8} = t_r4, {1,9} = t_b,
% {1,10} = r1_r1, {1.11} = r1_g1, {1,12} = r1_r2, {1,13} = r1_g2, 
% {1,14} = r1_r3, {1,15} = r1_g3, {1,16} = r1_r4,
% {1,17} = g1_g1, {1,18} = g1_r2, {1,19} = g1_g2, {1,20} = g1_r3,
% {1,21} = g1_g3, 
% {1,22} = r2_r2, {1,23} = r2_g2, {1,24} = r2_r3, {1,25} = g2_g2,
% {1,26} = r3_r3, {1,27} = g3_g3, {1,28} = r4_r4.

a = Geometry.a;
L = Geometry.RL(1,1);
g_Side = Geometry.g(1,1);
% g_Center = 0.22 : 0.02 : 0.3;
g_Center = Geometry.gCenterValue;

RadiationTable = cell(1,28);
TableParameter = cell(1,3);
RadiationAreaFactorTable = cell(1,28);

TablegaRatio = g_Center/a;
[Meshka, Meshga] = meshgrid(kaForTable,TablegaRatio);

TableParameter{1,1} = Meshka;
TableParameter{1,2} = Meshga;
TableParameter{1,3} = TablegaRatio;

%% Load Ztt
RadiationTable{1,1}=zeros(length(g_Center),length(kaForTable));

Area = pi*a^2; % Radiation Area
DivideFactor = rhowater*sp*Area;
RadiationAreaFactorTable{1,1}(1:length(g_Center),1) = DivideFactor;

% File Open
fid = fopen('tu_t.txt');

% Delete HeadLine
for i=1:5
    buffer = fgetl(fid);
end

% Load Data File
dataCell = textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f');
data = [dataCell{1} dataCell{2} dataCell{3} dataCell{4} dataCell{5}...
    dataCell{6}];

clear dataCell

% SplitData
Freq = data(:,1);
ZttData = data(:,2:end)./DivideFactor;

%ÆÄÀÏ ´Ý°í
fclose(fid);

% Rearrange
RadiationTable{1,1} = ZttData.';

%% Load Ztr1
RadiationTable{1,2}=zeros(length(g_Center),length(kaForTable));

Area = pi*a^2; % Radiation Area
DivideFactor = rhowater*sp*Area;
RadiationAreaFactorTable{1,2}(1:length(g_Center),1) = DivideFactor;

% File Open
fid = fopen('tu_r1.txt');

% Delete HeadLine
for i=1:5
    buffer = fgetl(fid);
end

% Load Data File
dataCell = textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f');
data = [dataCell{1} dataCell{2} dataCell{3} dataCell{4} dataCell{5}...
    dataCell{6}];

clear dataCell

% SplitData
Freq = data(:,1);
Ztr1Data = data(:,2:end)./DivideFactor;

%ÆÄÀÏ ´Ý°í
fclose(fid);

% Rearrange
RadiationTable{1,2} = Ztr1Data.';

%% Load Ztg1
RadiationTable{1,3}=zeros(length(g_Center),length(kaForTable));

Area = pi*a^2; % Radiation Area
DivideFactor = rhowater*sp*Area;
RadiationAreaFactorTable{1,3}(1:length(g_Center),1) = DivideFactor;

% File Open
fid = fopen('tu_g1.txt');

% Delete HeadLine
for i=1:5
    buffer = fgetl(fid);
end

% Load Data File
dataCell = textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f');
data = [dataCell{1} dataCell{2} dataCell{3} dataCell{4} dataCell{5}...
    dataCell{6}];

clear dataCell

% SplitData
Freq = data(:,1);
Ztg1Data = data(:,2:end)./DivideFactor;

%ÆÄÀÏ ´Ý°í
fclose(fid);

% Rearrange
RadiationTable{1,3} = Ztg1Data.';

%% Load Ztr2
RadiationTable{1,4}=zeros(length(g_Center),length(kaForTable));

Area = pi*a^2; % Radiation Area
DivideFactor = rhowater*sp*Area;
RadiationAreaFactorTable{1,4}(1:length(g_Center),1) = DivideFactor;

% File Open
fid = fopen('tu_r2.txt');

% Delete HeadLine
for i=1:5
    buffer = fgetl(fid);
end

% Load Data File
dataCell = textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f');
data = [dataCell{1} dataCell{2} dataCell{3} dataCell{4} dataCell{5}...
    dataCell{6}];

clear dataCell

% SplitData
Freq = data(:,1);
Ztr2Data = data(:,2:end)./DivideFactor;

%ÆÄÀÏ ´Ý°í
fclose(fid);

% Rearrange
RadiationTable{1,4} = Ztr2Data.';

%% Load Ztg2
RadiationTable{1,5}=zeros(length(g_Center),length(kaForTable));

Area = pi*a^2; % Radiation Area
DivideFactor = rhowater*sp*Area;
RadiationAreaFactorTable{1,5}(1:length(g_Center),1) = DivideFactor;

% File Open
fid = fopen('tu_g2.txt');

% Delete HeadLine
for i=1:5
    buffer = fgetl(fid);
end

% Load Data File
dataCell = textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f');
data = [dataCell{1} dataCell{2} dataCell{3} dataCell{4} dataCell{5}...
    dataCell{6}];

clear dataCell

% SplitData
Freq = data(:,1);
Ztg2Data = data(:,2:end)./DivideFactor;

%ÆÄÀÏ ´Ý°í
fclose(fid);

% Rearrange
RadiationTable{1,5} = Ztg2Data.';

%% Load Ztr3
RadiationTable{1,6}=zeros(length(g_Center),length(kaForTable));

Area = pi*a^2; % Radiation Area
DivideFactor = rhowater*sp*Area;
RadiationAreaFactorTable{1,6}(1:length(g_Center),1) = DivideFactor;

% File Open
fid = fopen('tu_r3.txt');

% Delete HeadLine
for i=1:5
    buffer = fgetl(fid);
end

% Load Data File
dataCell = textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f');
data = [dataCell{1} dataCell{2} dataCell{3} dataCell{4} dataCell{5}...
    dataCell{6}];

clear dataCell

% SplitData
Freq = data(:,1);
Ztr3Data = data(:,2:end)./DivideFactor;

%ÆÄÀÏ ´Ý°í
fclose(fid);

% Rearrange
RadiationTable{1,6} = Ztr3Data.';

%% Load Ztg3
RadiationTable{1,7}=zeros(length(g_Center),length(kaForTable));

Area = pi*a^2; % Radiation Area
DivideFactor = rhowater*sp*Area;
RadiationAreaFactorTable{1,7}(1:length(g_Center),1) = DivideFactor;

% File Open
fid = fopen('tu_g3.txt');

% Delete HeadLine
for i=1:5
    buffer = fgetl(fid);
end

% Load Data File
dataCell = textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f');
data = [dataCell{1} dataCell{2} dataCell{3} dataCell{4} dataCell{5}...
    dataCell{6}];

clear dataCell

% SplitData
Freq = data(:,1);
Ztg3Data = data(:,2:end)./DivideFactor;

%ÆÄÀÏ ´Ý°í
fclose(fid);

% Rearrange
RadiationTable{1,7} = Ztg3Data.';

%% Load Ztr4
RadiationTable{1,8}=zeros(length(g_Center),length(kaForTable));

Area = pi*a^2; % Radiation Area
DivideFactor = rhowater*sp*Area;
RadiationAreaFactorTable{1,8}(1:length(g_Center),1) = DivideFactor;

% File Open
fid = fopen('tu_r4.txt');

% Delete HeadLine
for i=1:5
    buffer = fgetl(fid);
end

% Load Data File
dataCell = textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f');
data = [dataCell{1} dataCell{2} dataCell{3} dataCell{4} dataCell{5}...
    dataCell{6}];

clear dataCell

% SplitData
Freq = data(:,1);
Ztr4Data = data(:,2:end)./DivideFactor;

%ÆÄÀÏ ´Ý°í
fclose(fid);

% Rearrange
RadiationTable{1,8} = Ztr4Data.';

%% Load Ztb
RadiationTable{1,9}=zeros(length(g_Center),length(kaForTable));

Area = pi*a^2; % Radiation Area
DivideFactor = rhowater*sp*Area;
RadiationAreaFactorTable{1,9}(1:length(g_Center),1) = DivideFactor;

% File Open
fid = fopen('tu_b.txt');

% Delete HeadLine
for i=1:5
    buffer = fgetl(fid);
end

% Load Data File
dataCell = textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f');
data = [dataCell{1} dataCell{2} dataCell{3} dataCell{4} dataCell{5}...
    dataCell{6}];

clear dataCell

% SplitData
Freq = data(:,1);
ZtbData = data(:,2:end)./DivideFactor;

%ÆÄÀÏ ´Ý°í
fclose(fid);

% Rearrange
RadiationTable{1,9} = ZtbData.';

%% Load Zr1r1
RadiationTable{1,10}=zeros(length(g_Center),length(kaForTable));

Area = 2*pi*a*L; % Radiation Area
DivideFactor = rhowater*sp*Area;
RadiationAreaFactorTable{1,10}(1:length(g_Center),1) = DivideFactor;

% File Open
fid = fopen('r1u_r1.txt');

% Delete HeadLine
for i=1:5
    buffer = fgetl(fid);
end

% Load Data File
dataCell = textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f');
data = [dataCell{1} dataCell{2} dataCell{3} dataCell{4} dataCell{5}...
    dataCell{6}];

clear dataCell

% SplitData
Freq = data(:,1);
Zr1r1Data = data(:,2:end)./DivideFactor;

%ÆÄÀÏ ´Ý°í
fclose(fid);

% Rearrange
RadiationTable{1,10} = Zr1r1Data.';

%% Load Zr1g1
RadiationTable{1,11}=zeros(length(g_Center),length(kaForTable));

Area = 2*pi*a*L; % Radiation Area
DivideFactor = rhowater*sp*Area;
RadiationAreaFactorTable{1,11}(1:length(g_Center),1) = DivideFactor;

% File Open
fid = fopen('r1u_g1.txt');

% Delete HeadLine
for i=1:5
    buffer = fgetl(fid);
end

% Load Data File
dataCell = textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f');
data = [dataCell{1} dataCell{2} dataCell{3} dataCell{4} dataCell{5}...
    dataCell{6}];

clear dataCell

% SplitData
Freq = data(:,1);
Zr1g1Data = data(:,2:end)./DivideFactor;

%ÆÄÀÏ ´Ý°í
fclose(fid);

% Rearrange
RadiationTable{1,11} = Zr1g1Data.';

%% Load Zr1r2
RadiationTable{1,12}=zeros(length(g_Center),length(kaForTable));

Area = 2*pi*a*L; % Radiation Area
DivideFactor = rhowater*sp*Area;
RadiationAreaFactorTable{1,12}(1:length(g_Center),1) = DivideFactor;

% File Open
fid = fopen('r1u_r2.txt');

% Delete HeadLine
for i=1:5
    buffer = fgetl(fid);
end

% Load Data File
dataCell = textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f');
data = [dataCell{1} dataCell{2} dataCell{3} dataCell{4} dataCell{5}...
    dataCell{6}];

clear dataCell

% SplitData
Freq = data(:,1);
Zr1r2Data = data(:,2:end)./DivideFactor;

%ÆÄÀÏ ´Ý°í
fclose(fid);

% Rearrange
RadiationTable{1,12} = Zr1r2Data.';

%% Load Zr1g2
RadiationTable{1,13}=zeros(length(g_Center),length(kaForTable));

Area = 2*pi*a*L; % Radiation Area
DivideFactor = rhowater*sp*Area;
RadiationAreaFactorTable{1,13}(1:length(g_Center),1) = DivideFactor;

% File Open
fid = fopen('r1u_g2.txt');

% Delete HeadLine
for i=1:5
    buffer = fgetl(fid);
end

% Load Data File
dataCell = textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f');
data = [dataCell{1} dataCell{2} dataCell{3} dataCell{4} dataCell{5}...
    dataCell{6}];

clear dataCell

% SplitData
Freq = data(:,1);
Zr1g2Data = data(:,2:end)./DivideFactor;

%ÆÄÀÏ ´Ý°í
fclose(fid);

% Rearrange
RadiationTable{1,13} = Zr1g2Data.';

%% Load Zr1r3
RadiationTable{1,14}=zeros(length(g_Center),length(kaForTable));

Area = 2*pi*a*L; % Radiation Area
DivideFactor = rhowater*sp*Area;
RadiationAreaFactorTable{1,14}(1:length(g_Center),1) = DivideFactor;

% File Open
fid = fopen('r1u_r3.txt');

% Delete HeadLine
for i=1:5
    buffer = fgetl(fid);
end

% Load Data File
dataCell = textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f');
data = [dataCell{1} dataCell{2} dataCell{3} dataCell{4} dataCell{5}...
    dataCell{6}];

clear dataCell

% SplitData
Freq = data(:,1);
Zr1r3Data = data(:,2:end)./DivideFactor;

%ÆÄÀÏ ´Ý°í
fclose(fid);

% Rearrange
RadiationTable{1,14} = Zr1r3Data.';

%% Load Zr1g3
RadiationTable{1,15}=zeros(length(g_Center),length(kaForTable));

Area = 2*pi*a*L; % Radiation Area
DivideFactor = rhowater*sp*Area;
RadiationAreaFactorTable{1,15}(1:length(g_Center),1) = DivideFactor;

% File Open
fid = fopen('r1u_g3.txt');

% Delete HeadLine
for i=1:5
    buffer = fgetl(fid);
end

% Load Data File
dataCell = textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f');
data = [dataCell{1} dataCell{2} dataCell{3} dataCell{4} dataCell{5}...
    dataCell{6}];

clear dataCell

% SplitData
Freq = data(:,1);
Zr1g3Data = data(:,2:end)./DivideFactor;

%ÆÄÀÏ ´Ý°í
fclose(fid);

% Rearrange
RadiationTable{1,15} = Zr1g3Data.';

%% Load Zr1r4
RadiationTable{1,16}=zeros(length(g_Center),length(kaForTable));

Area = 2*pi*a*L; % Radiation Area
DivideFactor = rhowater*sp*Area;
RadiationAreaFactorTable{1,16}(1:length(g_Center),1) = DivideFactor;

% File Open
fid = fopen('r1u_r4.txt');

% Delete HeadLine
for i=1:5
    buffer = fgetl(fid);
end

% Load Data File
dataCell = textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f');
data = [dataCell{1} dataCell{2} dataCell{3} dataCell{4} dataCell{5}...
    dataCell{6}];

clear dataCell

% SplitData
Freq = data(:,1);
Zr1r4Data = data(:,2:end)./DivideFactor;

%ÆÄÀÏ ´Ý°í
fclose(fid);

% Rearrange
RadiationTable{1,16} = Zr1r4Data.';

%% Load Zg1g1
RadiationTable{1,17}=zeros(length(g_Center),length(kaForTable));

Area = 2*pi*a*g_Side; % Radiation Area
DivideFactor = rhowater*sp*Area;
RadiationAreaFactorTable{1,17}(1:length(g_Center),1) = DivideFactor;

% File Open
fid = fopen('g1u_g1.txt');

% Delete HeadLine
for i=1:5
    buffer = fgetl(fid);
end

% Load Data File
dataCell = textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f');
data = [dataCell{1} dataCell{2} dataCell{3} dataCell{4} dataCell{5}...
    dataCell{6}];

clear dataCell

% SplitData
Freq = data(:,1);
Zg1g1Data = data(:,2:end)./DivideFactor;

%ÆÄÀÏ ´Ý°í
fclose(fid);

% Rearrange
RadiationTable{1,17} = Zg1g1Data.';

%% Load Zg1r2
RadiationTable{1,18}=zeros(length(g_Center),length(kaForTable));

Area = 2*pi*a*g_Side; % Radiation Area
DivideFactor = rhowater*sp*Area;
RadiationAreaFactorTable{1,18}(1:length(g_Center),1) = DivideFactor;

% File Open
fid = fopen('g1u_r2.txt');

% Delete HeadLine
for i=1:5
    buffer = fgetl(fid);
end

% Load Data File
dataCell = textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f');
data = [dataCell{1} dataCell{2} dataCell{3} dataCell{4} dataCell{5}...
    dataCell{6}];

clear dataCell

% SplitData
Freq = data(:,1);
Zg1r2Data = data(:,2:end)./DivideFactor;

%ÆÄÀÏ ´Ý°í
fclose(fid);

% Rearrange
RadiationTable{1,18} = Zg1r2Data.';

%% Load Zg1g2
RadiationTable{1,19}=zeros(length(g_Center),length(kaForTable));

Area = 2*pi*a*g_Side; % Radiation Area
DivideFactor = rhowater*sp*Area;
RadiationAreaFactorTable{1,19}(1:length(g_Center),1) = DivideFactor;

% File Open
fid = fopen('g1u_g2.txt');

% Delete HeadLine
for i=1:5
    buffer = fgetl(fid);
end

% Load Data File
dataCell = textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f');
data = [dataCell{1} dataCell{2} dataCell{3} dataCell{4} dataCell{5}...
    dataCell{6}];

clear dataCell

% SplitData
Freq = data(:,1);
Zg1g2Data = data(:,2:end)./DivideFactor;

%ÆÄÀÏ ´Ý°í
fclose(fid);

% Rearrange
RadiationTable{1,19} = Zg1g2Data.';

%% Load Zg1r3
RadiationTable{1,20}=zeros(length(g_Center),length(kaForTable));

Area = 2*pi*a*g_Side; % Radiation Area
DivideFactor = rhowater*sp*Area;
RadiationAreaFactorTable{1,20}(1:length(g_Center),1) = DivideFactor;

% File Open
fid = fopen('g1u_r3.txt');

% Delete HeadLine
for i=1:5
    buffer = fgetl(fid);
end

% Load Data File
dataCell = textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f');
data = [dataCell{1} dataCell{2} dataCell{3} dataCell{4} dataCell{5}...
    dataCell{6}];

clear dataCell

% SplitData
Freq = data(:,1);
Zg1r3Data = data(:,2:end)./DivideFactor;

%ÆÄÀÏ ´Ý°í
fclose(fid);

% Rearrange
RadiationTable{1,20} = Zg1r3Data.';

%% Load Zg1g3
RadiationTable{1,21}=zeros(length(g_Center),length(kaForTable));

Area = 2*pi*a*g_Side; % Radiation Area
DivideFactor = rhowater*sp*Area;
RadiationAreaFactorTable{1,21}(1:length(g_Center),1) = DivideFactor;

% File Open
fid = fopen('g1u_g3.txt');

% Delete HeadLine
for i=1:5
    buffer = fgetl(fid);
end

% Load Data File
dataCell = textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f');
data = [dataCell{1} dataCell{2} dataCell{3} dataCell{4} dataCell{5}...
    dataCell{6}];

clear dataCell

% SplitData
Freq = data(:,1);
Zg1g3Data = data(:,2:end)./DivideFactor;

%ÆÄÀÏ ´Ý°í
fclose(fid);

% Rearrange
RadiationTable{1,21} = Zg1g3Data.';

%% Load Zr2r2
RadiationTable{1,22}=zeros(length(g_Center),length(kaForTable));

Area = 2*pi*a*L; % Radiation Area
DivideFactor = rhowater*sp*Area;
RadiationAreaFactorTable{1,22}(1:length(g_Center),1) = DivideFactor;

% File Open
fid = fopen('r2u_r2.txt');

% Delete HeadLine
for i=1:5
    buffer = fgetl(fid);
end

% Load Data File
dataCell = textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f');
data = [dataCell{1} dataCell{2} dataCell{3} dataCell{4} dataCell{5}...
    dataCell{6}];

clear dataCell

% SplitData
Freq = data(:,1);
Zr2r2Data = data(:,2:end)./DivideFactor;

%ÆÄÀÏ ´Ý°í
fclose(fid);

% Rearrange
RadiationTable{1,22} = Zr2r2Data.';

%% Load Zr2g2
RadiationTable{1,23}=zeros(length(g_Center),length(kaForTable));

Area = 2*pi*a*L; % Radiation Area
DivideFactor = rhowater*sp*Area;
RadiationAreaFactorTable{1,23}(1:length(g_Center),1) = DivideFactor;

% File Open
fid = fopen('r2u_g2.txt');

% Delete HeadLine
for i=1:5
    buffer = fgetl(fid);
end

% Load Data File
dataCell = textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f');
data = [dataCell{1} dataCell{2} dataCell{3} dataCell{4} dataCell{5}...
    dataCell{6}];

clear dataCell

% SplitData
Freq = data(:,1);
Zr2g2Data = data(:,2:end)./DivideFactor;

%ÆÄÀÏ ´Ý°í
fclose(fid);

% Rearrange
RadiationTable{1,23} = Zr2g2Data.';

%% Load Zr2r3
RadiationTable{1,24}=zeros(length(g_Center),length(kaForTable));

Area = 2*pi*a*L; % Radiation Area
DivideFactor = rhowater*sp*Area;
RadiationAreaFactorTable{1,24}(1:length(g_Center),1) = DivideFactor;

% File Open
fid = fopen('r2u_r3.txt');

% Delete HeadLine
for i=1:5
    buffer = fgetl(fid);
end

% Load Data File
dataCell = textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f');
data = [dataCell{1} dataCell{2} dataCell{3} dataCell{4} dataCell{5}...
    dataCell{6}];

clear dataCell

% SplitData
Freq = data(:,1);
Zr2r3Data = data(:,2:end)./DivideFactor;

%ÆÄÀÏ ´Ý°í
fclose(fid);

% Rearrange
RadiationTable{1,24} = Zr2r3Data.';

%% Load Zg2g2
RadiationTable{1,25}=zeros(length(g_Center),length(kaForTable));

Area = 2*pi*a*g_Center; % Radiation Area
DivideFactor = rhowater*sp*Area;
RadiationAreaFactorTable{1,25}(1:length(g_Center),1) = DivideFactor;

% File Open
fid = fopen('g2u_g2.txt');

% Delete HeadLine
for i=1:5
    buffer = fgetl(fid);
end

% Load Data File
dataCell = textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f');
data = [dataCell{1} dataCell{2} dataCell{3} dataCell{4} dataCell{5}...
    dataCell{6}];

clear dataCell

% SplitData
Freq = data(:,1);
for ii = 1 : length(g_Center);
    Zg2g2Data(:,ii) = data(:,ii+1)./DivideFactor(ii);
end

%ÆÄÀÏ ´Ý°í
fclose(fid);

% Rearrange
RadiationTable{1,25} = Zg2g2Data.';

%% Load Zr3r3
RadiationTable{1,26}=zeros(length(g_Center),length(kaForTable));

Area = 2*pi*a*L; % Radiation Area
DivideFactor = rhowater*sp*Area;
RadiationAreaFactorTable{1,26}(1:length(g_Center),1) = DivideFactor;

% File Open
fid = fopen('r3u_r3.txt');

% Delete HeadLine
for i=1:5
    buffer = fgetl(fid);
end

% Load Data File
dataCell = textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f');
data = [dataCell{1} dataCell{2} dataCell{3} dataCell{4} dataCell{5}...
    dataCell{6}];

clear dataCell

% SplitData
Freq = data(:,1);
Zr3r3Data = data(:,2:end)./DivideFactor;

%ÆÄÀÏ ´Ý°í
fclose(fid);

% Rearrange
RadiationTable{1,26} = Zr3r3Data.';

%% Load Zg3g3
RadiationTable{1,27}=zeros(length(g_Center),length(kaForTable));

Area = 2*pi*a*g_Side; % Radiation Area
DivideFactor = rhowater*sp*Area;
RadiationAreaFactorTable{1,27}(1:length(g_Center),1) = DivideFactor;

% File Open
fid = fopen('g3u_g3.txt');

% Delete HeadLine
for i=1:5
    buffer = fgetl(fid);
end

% Load Data File
dataCell = textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f');
data = [dataCell{1} dataCell{2} dataCell{3} dataCell{4} dataCell{5}...
    dataCell{6}];

clear dataCell

% SplitData
Freq = data(:,1);
Zg3g3Data = data(:,2:end)./DivideFactor;

%ÆÄÀÏ ´Ý°í
fclose(fid);

% Rearrange
RadiationTable{1,27} = Zg3g3Data.';

%% Load Zr4r4
RadiationTable{1,28}=zeros(length(g_Center),length(kaForTable));

Area = 2*pi*a*L; % Radiation Area
DivideFactor = rhowater*sp*Area;
RadiationAreaFactorTable{1,28}(1:length(g_Center),1) = DivideFactor;

% File Open
fid = fopen('r4u_r4.txt');

% Delete HeadLine
for i=1:5
    buffer = fgetl(fid);
end

% Load Data File
dataCell = textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f');
data = [dataCell{1} dataCell{2} dataCell{3} dataCell{4} dataCell{5}...
    dataCell{6}];

clear dataCell

% SplitData
Freq = data(:,1);
Zr4r4Data = data(:,2:end)./DivideFactor;

%ÆÄÀÏ ´Ý°í
fclose(fid);

% Rearrange
RadiationTable{1,28} = Zr4r4Data.';

%% Load Zgg Sample Driving g
% RadiationTable{1,7}=zeros(length(g_side),length(kaForTable));
% 
% Area = 2*pi*a*g_side; % Radiation Area
% DivideFactor = rhowater*sp*Area;
% RadiationAreaFactorTable{1,7}(1:length(g_side),1) = DivideFactor;
% 
% % File Open
% fid = fopen('gu_g.txt');
% 
% % Delete HeadLine
% for i=1:5
%     buffer = fgetl(fid);
% end
% 
% % Load Data File
% dataCell = textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f');
% data = [dataCell{1} dataCell{2} dataCell{3} dataCell{4} dataCell{5}...
%     dataCell{6} dataCell{7} dataCell{8} dataCell{9} dataCell{10} dataCell{11}];
% 
% clear dataCell
% 
% % SplitData
% Freq = data(:,1);
% for ii = 1 : length(g_side);
%     ZggData(:,ii) = data(:,ii+1)./DivideFactor(ii);
% end
% 
% %ÆÄÀÏ ´Ý°í
% fclose(fid);
% 
% % Rearrange
% RadiationTable{1,7} = ZggData.';
end