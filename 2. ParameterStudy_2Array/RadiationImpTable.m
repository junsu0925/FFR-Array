function [RadiationTable,RadiationAreaFactorTable,TableParameter] = RadiationImpTable(a, L, gmpnum, rhowater, sp, kaForTable)
% {1,1} = r1_b, {1.2} = r1_r1, {1,3} = r1_g, {1,4} = r1_r2, {1,5} = r1_t
% {1,6} = g_b {1,7} = g_g1 {1,8} = b_b {1,9} = b_t

RadiationTable = cell(1,9);
TableParameter = cell(1,3);
RadiationAreaFactorTable = cell(1,9);

TablegaRatio = gmpnum/a;
[Meshka, Meshga] = meshgrid(kaForTable,TablegaRatio);

TableParameter{1,1} = Meshka;
TableParameter{1,2} = Meshga;
TableParameter{1,3} = TablegaRatio;

%% Load Zr1b
RadiationTable{1,1}=zeros(length(gmpnum),length(kaForTable));

Area = 2*pi*a*L; % Radiation Area
DivideFactor = rhowater*sp*Area;
RadiationAreaFactorTable{1,1}(1:length(gmpnum),1) = DivideFactor;

% File Open
fid = fopen('r1u_b.txt');

% Delete HeadLine
for i=1:5
    buffer = fgetl(fid);
end

% Load Data File
dataCell = textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f');
data = [dataCell{1} dataCell{2} dataCell{3} dataCell{4} dataCell{5}...
    dataCell{6} dataCell{7} dataCell{8} dataCell{9} dataCell{10} dataCell{11}];

clear dataCell

% SplitData
Freq = data(:,1);
Zr1bData = data(:,2:end)./DivideFactor;

%ÆÄÀÏ ´Ý°í
fclose(fid);

% Rearrange
RadiationTable{1,1} = Zr1bData.';

%% Load Zr1r1
RadiationTable{1,2}=zeros(length(gmpnum),length(kaForTable));

Area = 2*pi*a*L; % Radiation Area
DivideFactor = rhowater*sp*Area;
RadiationAreaFactorTable{1,2}(1:length(gmpnum),1) = DivideFactor;

% File Open
fid = fopen('r1u_r1.txt');

% Delete HeadLine
for i=1:5
    buffer = fgetl(fid);
end

% Load Data File
dataCell = textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f');
data = [dataCell{1} dataCell{2} dataCell{3} dataCell{4} dataCell{5}...
    dataCell{6} dataCell{7} dataCell{8} dataCell{9} dataCell{10} dataCell{11}];

clear dataCell

% SplitData
Freq = data(:,1);
Zr1r1Data = data(:,2:end)./DivideFactor;

%ÆÄÀÏ ´Ý°í
fclose(fid);

% Rearrange
RadiationTable{1,2} = Zr1r1Data.';

%% Load Zr1g
RadiationTable{1,3}=zeros(length(gmpnum),length(kaForTable));

Area = 2*pi*a*L; % Radiation Area
DivideFactor = rhowater*sp*Area;
RadiationAreaFactorTable{1,3}(1:length(gmpnum),1) = DivideFactor;

% File Open
fid = fopen('r1u_g.txt');

% Delete HeadLine
for i=1:5
    buffer = fgetl(fid);
end

% Load Data File
dataCell = textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f');
data = [dataCell{1} dataCell{2} dataCell{3} dataCell{4} dataCell{5}...
    dataCell{6} dataCell{7} dataCell{8} dataCell{9} dataCell{10} dataCell{11}];

clear dataCell

% SplitData
Freq = data(:,1);
Zr1gData = data(:,2:end)./DivideFactor;

%ÆÄÀÏ ´Ý°í
fclose(fid);

% Rearrange
RadiationTable{1,3} = Zr1gData.';

%% Load Zr1r2
RadiationTable{1,4}=zeros(length(gmpnum),length(kaForTable));

Area = 2*pi*a*L; % Radiation Area
DivideFactor = rhowater*sp*Area;
RadiationAreaFactorTable{1,4}(1:length(gmpnum),1) = DivideFactor;

% File Open
fid = fopen('r1u_r2.txt');

% Delete HeadLine
for i=1:5
    buffer = fgetl(fid);
end

% Load Data File
dataCell = textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f');
data = [dataCell{1} dataCell{2} dataCell{3} dataCell{4} dataCell{5}...
    dataCell{6} dataCell{7} dataCell{8} dataCell{9} dataCell{10} dataCell{11}];

clear dataCell

% SplitData
Freq = data(:,1);
Zr1r2Data = data(:,2:end)./DivideFactor;

%ÆÄÀÏ ´Ý°í
fclose(fid);

% Rearrange
RadiationTable{1,4} = Zr1r2Data.';

%% Load Zr1t
RadiationTable{1,5}=zeros(length(gmpnum),length(kaForTable));

Area = 2*pi*a*L; % Radiation Area
DivideFactor = rhowater*sp*Area;
RadiationAreaFactorTable{1,5}(1:length(gmpnum),1) = DivideFactor;

% File Open
fid = fopen('r1u_t.txt');

% Delete HeadLine
for i=1:5
    buffer = fgetl(fid);
end

% Load Data File
dataCell = textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f');
data = [dataCell{1} dataCell{2} dataCell{3} dataCell{4} dataCell{5}...
    dataCell{6} dataCell{7} dataCell{8} dataCell{9} dataCell{10} dataCell{11}];

clear dataCell

% SplitData
Freq = data(:,1);
Zr1tData = data(:,2:end)./DivideFactor;

%ÆÄÀÏ ´Ý°í
fclose(fid);

% Rearrange
RadiationTable{1,5} = Zr1tData.';

%% Load Zgb
RadiationTable{1,6}=zeros(length(gmpnum),length(kaForTable));

Area = 2*pi*a*gmpnum; % Radiation Area
DivideFactor = rhowater*sp*Area;
RadiationAreaFactorTable{1,6}(1:length(gmpnum),1) = DivideFactor;

% File Open
fid = fopen('gu_b.txt');

% Delete HeadLine
for i=1:5
    buffer = fgetl(fid);
end

% Load Data File
dataCell = textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f');
data = [dataCell{1} dataCell{2} dataCell{3} dataCell{4} dataCell{5}...
    dataCell{6} dataCell{7} dataCell{8} dataCell{9} dataCell{10} dataCell{11}];

clear dataCell

% SplitData
Freq = data(:,1);
for ii = 1 : length(gmpnum);
    ZgbData(:,ii) = data(:,ii+1)./DivideFactor(ii);
end

%ÆÄÀÏ ´Ý°í
fclose(fid);

% Rearrange
RadiationTable{1,6} = ZgbData.';

%% Load Zgg
RadiationTable{1,7}=zeros(length(gmpnum),length(kaForTable));

Area = 2*pi*a*gmpnum; % Radiation Area
DivideFactor = rhowater*sp*Area;
RadiationAreaFactorTable{1,7}(1:length(gmpnum),1) = DivideFactor;

% File Open
fid = fopen('gu_g.txt');

% Delete HeadLine
for i=1:5
    buffer = fgetl(fid);
end

% Load Data File
dataCell = textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f');
data = [dataCell{1} dataCell{2} dataCell{3} dataCell{4} dataCell{5}...
    dataCell{6} dataCell{7} dataCell{8} dataCell{9} dataCell{10} dataCell{11}];

clear dataCell

% SplitData
Freq = data(:,1);
for ii = 1 : length(gmpnum);
    ZggData(:,ii) = data(:,ii+1)./DivideFactor(ii);
end

%ÆÄÀÏ ´Ý°í
fclose(fid);

% Rearrange
RadiationTable{1,7} = ZggData.';

%% Load Zbb
RadiationTable{1,8}=zeros(length(gmpnum),length(kaForTable));

Area = pi*a^2; % Radiation Area
DivideFactor = rhowater*sp*Area;
RadiationAreaFactorTable{1,8}(1:length(gmpnum),1) = DivideFactor;

% File Open
fid = fopen('bu_b.txt');

% Delete HeadLine
for i=1:5
    buffer = fgetl(fid);
end

% Load Data File
dataCell = textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f');
data = [dataCell{1} dataCell{2} dataCell{3} dataCell{4} dataCell{5}...
    dataCell{6} dataCell{7} dataCell{8} dataCell{9} dataCell{10} dataCell{11}];

clear dataCell

% SplitData
Freq = data(:,1);
ZbbData = data(:,2:end)./DivideFactor;

%ÆÄÀÏ ´Ý°í
fclose(fid);

% Rearrange
RadiationTable{1,8} = ZbbData.';

%% Load Zbt
RadiationTable{1,9}=zeros(length(gmpnum),length(kaForTable));

Area = pi*a^2; % Radiation Area
DivideFactor = rhowater*sp*Area;
RadiationAreaFactorTable{1,9}(1:length(gmpnum),1) = DivideFactor;

% File Open
fid = fopen('bu_t.txt');

% Delete HeadLine
for i=1:5
    buffer = fgetl(fid);
end

% Load Data File
dataCell = textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f');
data = [dataCell{1} dataCell{2} dataCell{3} dataCell{4} dataCell{5}...
    dataCell{6} dataCell{7} dataCell{8} dataCell{9} dataCell{10} dataCell{11}];

clear dataCell

% SplitData
Freq = data(:,1);
ZbtData = data(:,2:end)./DivideFactor;

%ÆÄÀÏ ´Ý°í
fclose(fid);

% Rearrange
RadiationTable{1,9} = ZbtData.';

end