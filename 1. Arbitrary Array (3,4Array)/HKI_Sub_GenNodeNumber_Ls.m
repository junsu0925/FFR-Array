function [NN] = HKI_Sub_GenNodeNumber_Ls(MagDVector,freq,c0)
% Freq �������� Node ���� ���
SliceNum = 21;
MinNodeNumber = 21;

RefCalLen = (c0/freq)/SliceNum;
NN_temp1 = ceil(MagDVector./RefCalLen);
OddEven = find(mod(NN_temp1,2)==0);
MakeOddNumber = zeros(1,length(NN_temp1));
MakeOddNumber(OddEven) = 1;

NN = NN_temp1+MakeOddNumber;

% Ư�� ������ ����
% NN = MinNodeNumber*ones(1,length(NN));

% Node ������ 5�� �̸��� ���, 5���� ����
if min(NN) < MinNodeNumber
    SliceNum = (c0/freq)/(min(MagDVector)/MinNodeNumber);
    RefCalLen = (c0/freq)/SliceNum;
    NN_temp1 = ceil(MagDVector./RefCalLen);
    OddEven = find(mod(NN_temp1,2)==0);
    MakeOddNumber = zeros(1,length(NN_temp1));
    MakeOddNumber(OddEven) = 1;
    
    NN = NN_temp1+MakeOddNumber;
end

% if NN(1) < MinNodeNumber
%     SliceNum = (c0/freq)/(MagDVector(1)/MinNodeNumber);
%     RefCalLen = (c0/freq)/SliceNum;
%     NN_temp1 = ceil(MagDVector./RefCalLen);
%     OddEven = find(mod(NN_temp1,2)==0);
%     MakeOddNumber = zeros(1,length(NN_temp1));
%     MakeOddNumber(OddEven) = 1;
%     
%     NN = NN_temp1+MakeOddNumber;
% end

for VuNum = 1:length(NN)
    SFName = sprintf('Surface %d', VuNum);
    SurfaceNumber(VuNum,1) = {SFName};
end
NumberOfNodes = NN';
disp('The node number assigned to the each surface is as follows.')
disp(' ')
FreqDisp = sprintf('         Frequency =  %d [Hz]', freq);
disp(FreqDisp)
tblA = table(SurfaceNumber,NumberOfNodes);
disp(tblA)
end