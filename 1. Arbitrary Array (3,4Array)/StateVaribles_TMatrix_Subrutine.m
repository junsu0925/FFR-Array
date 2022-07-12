function [StateVector, TMatrix]=StateVaribles_TMatrix_Subrutine(zlength1,zlength2,N)
%% State variables of the system
StateVector.Ua = zeros(zlength1,1); % Node Velocity vector
StateVector.Fa= zeros(zlength1,1); % Node Force vector
StateVector.RVa = zeros(N,1); % Ring Voltage vector
StateVector.RCa = zeros(N,1); % Ring Current vector

StateVector.Ui = zeros(3,1,N); % Node Velocity vector for i-th Ring
StateVector.Fi = zeros(3,1,N); % Node Force vector for i-th Ring
StateVector.RVi = zeros(1,1,N); % Applied Voltage vector for i-th Ring
StateVector.RCi = zeros(1,1,N); % Applied Current vector for i-th Ring

StateVector.Ug = zeros(3,1,N-1); % Gab Velocity vector for i-th Ring
StateVector.Fg = zeros(3,1,N-1); % Gab Force vector for i-th Ring

StateVector.Uas = zeros(zlength2,1); % Acoustic Raidation Surface Velocity vector

%% Transformation Matrix of the system
% Transformation Matrix from the velocity vector of the global system of
% the PZT Ring
TMatrix.Tsre = zeros(3,zlength1,N); 
for Count1 = 1 : N % i
    for Count2 = 1 : 3 % n
        for Count3 = 1 : zlength1 % m
            if (Count3 == 4 * Count1 - 4 + Count2)
                TMatrix.Tsre(Count2,Count3,Count1) = 1;
            end
        end
    end
end
clear Count1 Count2 Count3

% Transformation Matrix from the velocity vector of the global system of
% the Gab
TMatrix.Tsge = zeros(3,zlength1,N-1); 
for Count1 = 1 : N-1 % i
    for Count2 = 1 : 3 % n
        for Count3 = 1 : zlength1 % m
            if (Count3 == 4 * Count1 - 2 + Count2)
                TMatrix.Tsge(Count2,Count3,Count1) = 1;
            end
        end
    end
end
clear Count1 Count2 Count3

% Transformation Matrix from the velocity vector of the global system of
% acoustic radiation surface
TMatrix.Tsa = zeros(zlength2,zlength1); 
for Count1 = 1 : zlength2 % n
    for Count2 = 1 : zlength1 % m
        if Count1 == 1 && Count2 == 1
            TMatrix.Tsa(Count1,Count2) = 1;
        elseif (Count1 == 2*N+1) && (Count2 == zlength1)
            TMatrix.Tsa(Count1,Count2) = 1;
        elseif (Count2 == 2*(Count1-1)) && (2 <= Count1 <= 2*N)
            TMatrix.Tsa(Count1,Count2) = 1;
        end
    end
end
clear Count1 Count2

% Transformation Matrix from the velocity vector of the global system of
% the elevtrical state vactors of PZT ring
TMatrix.Tc = zeros(zlength1,N);
for Count1 = 1 : zlength1 % m
    for Count2 = 1 : N % n
        if (Count1 == 4*Count2-2)
            TMatrix.Tc(Count1,Count2) = 1;
        end
    end
end
clear Count1 Count2
end


% %% State variables of the system
% StateVector.Ua = zeros(zlength1,1); % Node Velocity vector
% StateVector.Fa= zeros(zlength1,1); % Node Force vector
% StateVector.RVa = zeros(N,1); % Ring Voltage vector
% StateVector.RCa = zeros(N,1); % Ring Current vector
% 
% StateVector.Ui = zeros(3,1,N); % Node Velocity vector for i-th Ring
% StateVector.Fi = zeros(3,1,N); % Node Force vector for i-th Ring
% StateVector.RVi = zeros(1,1,N); % Applied Voltage vector for i-th Ring
% StateVector.RCi = zeros(1,1,N); % Applied Current vector for i-th Ring
% 
% StateVector.Ug = zeros(3,1,N-1); % Gab Velocity vector for i-th Ring
% StateVector.Fg = zeros(3,1,N-1); % Gab Force vector for i-th Ring
% 
% StateVector.Uas = zeros(zlength2,1); % Acoustic Raidation Surface Velocity vector
% 
% %% Transformation Matrix of the system
% % Transformation Matrix from the velocity vector of the global system of
% % the PZT Ring
% TMatrix.TRsre = zeros(3,zlength1,N); 
% for Count1 = 1 : N % i
%     for Count2 = 1 : 3 % n
%         for Count3 = 1 : zlength1 % m
%             if (Count3 == 4 * Count1 - 4 + Count2)
%                 TMatrix.TRsre(Count2,Count3,Count1) = 1;
%             end
%         end
%     end
% end
% clear Count1 Count2 Count3
% 
% % Transformation Matrix from the velocity vector of the global system of
% % the Gab
% TMatrix.Tgsre = zeros(3,zlength1,N-1); 
% for Count1 = 1 : N-1 % i
%     for Count2 = 1 : 3 % n
%         for Count3 = 1 : zlength1 % m
%             if (Count3 == 4 * Count1 - 2 + Count2)
%                 TMatrix.Tgsre(Count2,Count3,Count1) = 1;
%             end
%         end
%     end
% end
% clear Count1 Count2 Count3
% 
% % Transformation Matrix from the velocity vector of the global system of
% % acoustic radiation surface
% TMatrix.Tsa = zeros(zlength2,zlength1); 
% for Count1 = 1 : zlength2 % n
%     for Count2 = 1 : zlength1 % m
%         if Count1 == 1 && Count2 == 1
%             TMatrix.Tsa(Count1,Count2) = 1;
%         elseif (Count1 == 2*N+1) && (Count2 == zlength1)
%             TMatrix.Tsa(Count1,Count2) = 1;
%         elseif (Count2 == 2*(Count1-1)) && (2 <= Count1 <= 2*N)
%             TMatrix.Tsa(Count1,Count2) = 1;
%         end
%     end
% end
% clear Count1 Count2
% 
% % Transformation Matrix from the velocity vector of the global system of
% % the elevtrical state vactors of PZT ring
% TMatrix.Tc = zeros(zlength1,N);
% for Count1 = 1 : zlength1 % m
%     for Count2 = 1 : N % n
%         if (Count1 == 4*Count2-2)
%             TMatrix.Tc(Count1,Count2) = 1;
%         end
%     end
% end