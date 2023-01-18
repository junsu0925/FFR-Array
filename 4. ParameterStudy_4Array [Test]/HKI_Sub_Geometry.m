function [Line_Sec,RealAcuteAngle,MagDVector] = HKI_Sub_Geometry(Geometry,Num)

a = Geometry.a; % Ring Average Radius
l = Geometry.RL(1);
g = Geometry.g(1);

if Num == 1
   Line_Sec(:,:,1) = [0 -l/2; a -l/2];
   Line_Sec(:,:,2) = [a -l/2; a l/2];
   Line_Sec(:,:,3) = [a l/2; 0 l/2];
elseif Num == 2
    Line_Sec(:,:,1) = [0 -l-g/2;a -l-g/2];
    Line_Sec(:,:,2) = [a -l-g/2;a -g/2];
    Line_Sec(:,:,3) = [a -g/2;a g/2];
    Line_Sec(:,:,4) = [a g/2;a l+g/2];
    Line_Sec(:,:,5) = [a l+g/2;0 l+g/2];
elseif Num == 3
    Line_Sec(:,:,1) = [0 -3/2*l-g;a -3/2*l-g];
    Line_Sec(:,:,2) = [a -3/2*l-g;a -1/2*l-g];
    Line_Sec(:,:,3) = [a -1/2*l-g;a -1/2*l];
    Line_Sec(:,:,4) = [a -1/2*l;a 1/2*l];
    Line_Sec(:,:,5) = [a 1/2*l;a 1/2*l+g];
    Line_Sec(:,:,6) = [a 1/2*l+g;a 3/2*l+g];
    Line_Sec(:,:,7) = [a 3/2*l+g;0 3/2*l+g];
elseif Num == 4
    g1 = Geometry.g(1);
    g2 = Geometry.g(2);
    g3 = Geometry.g(3);

    Line_Sec(:,:,1) = [0 -2*l-g1-1/2*g2;a -2*l-g1-1/2*g2];
    Line_Sec(:,:,2) = [a -2*l-g1-1/2*g2;a -l-g1-1/2*g2];
    Line_Sec(:,:,3) = [a -l-g1-1/2*g2;a -l-1/2*g2];
    Line_Sec(:,:,4) = [a -l-1/2*g2;a -1/2*g2];
    Line_Sec(:,:,5) = [a -1/2*g2;a 1/2*g2];
    Line_Sec(:,:,6) = [a 1/2*g2;a l+1/2*g2];
    Line_Sec(:,:,7) = [a l+1/2*g2;a l+1/2*g2+g3];
    Line_Sec(:,:,8) = [a l+1/2*g2+g3;a 2*l+1/2*g2+g3];
    Line_Sec(:,:,9) = [a 2*l+1/2*g2+g3;0 2*l+1/2*g2+g3];

%     Line_Sec(:,:,1) = [0 -2*l-3/2*g;a -2*l-3/2*g];
%     Line_Sec(:,:,2) = [a -2*l-3/2*g;a -l-3/2*g];
%     Line_Sec(:,:,3) = [a -l-3/2*g;a -l-1/2*g];
%     Line_Sec(:,:,4) = [a -l-1/2*g;a -1/2*g];
%     Line_Sec(:,:,5) = [a -1/2*g;a 1/2*g];
%     Line_Sec(:,:,6) = [a 1/2*g;a l+1/2*g];
%     Line_Sec(:,:,7) = [a l+1/2*g;a l+3/2*g];
%     Line_Sec(:,:,8) = [a l+3/2*g;a 2*l+3/2*g];
%     Line_Sec(:,:,9) = [a 2*l+3/2*g;0 2*l+3/2*g];
end

%% Line 정의 및 Normal Vector 검증
[VerifiedLine,RealAcuteAngle,MagDVector] = HKI_Sub_VerifiedLine_Ls(Line_Sec);

end
