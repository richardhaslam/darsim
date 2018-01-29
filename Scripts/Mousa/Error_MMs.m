close all; clear all; clc;
Directory = '../../Results_Comparison/Input/';

N_TestCase = 6;
TestCase_Name = 'SinglePhase';

Lx = 270;  Nx = 270;  dx = Lx/Nx;
Ly = 270;  Ny = 270;  dy = Ly/Ny;
Lz = 0.1;  Nz = 001;  dz = Lz/Nz;
 
xcm = linspace(dx/2 , Lx-dx/2 , Nx)';
ycm = linspace(dy/2 , Ly-dy/2 , Ny)';
zcm = linspace(dz/2 , Lz-dz/2 , Nz)';

Nt = 1;    % Number of time-steps
dt = 1e-4; % Size of time-step [day]

Error_P = zeros(N_TestCase-1,Nt);
Error_S = zeros(N_TestCase-1,Nt);

%% Looping over time
for t = 1:Nt
    P = zeros(Nx*Ny*Nz,N_TestCase);
    S = zeros(Nx*Ny*Nz,N_TestCase);
    
    fprintf('Checking timestep: %f\n', t);
    % Reading Files
    for n = 1:N_TestCase
        File = strcat(Directory, TestCase_Name,'_L',num2str(n-1),'_Sol',num2str(t),'.txt');
        FileID = fopen(File, 'r');
        Matrix = textscan(FileID, '%s');
        Matrix = str2double(Matrix{1,1});
        Matrix = reshape(Matrix,[3,Nx*Ny*Nz])';
        fclose(FileID);
        P(:,n) = Matrix(:,2);
        S(:,n) = Matrix(:,3);
    end
    for n = 1:N_TestCase - 1
        Error_P(n,t) = norm( (P(:,1) - P(:,n+1)) ) / length(P(:,1));
        Error_S(n,t) = norm( (S(:,1) - S(:,n+1)) ) / length(S(:,1));
    end
end

% %% Ploting Pressure Error
% figure(2);
% legendInfo = cell(length(Nx)-1,1);
% for n=1:length(Nx)-1
%     legendInfo{n} = strcat('DNS\_',num2str(Nx(1)),'\times',num2str(Ny(1)),' vs EDFM\_',num2str(Nx(n+1)),'\times',num2str(Ny(n+1)));
%     plot( (1:t)*dt , Error_P(:,n) ,'LineWidth',2 );
%     hold on;
% end
% title('Pressure Error');
% legend(legendInfo);
% xlabel('Simulation Time [day], (\Deltat = 10^{-4} day)');
% ylabel('Error(P)');
% 
% %% Ploting Saturation Error
% figure(3);
% legendInfo = cell(length(Nx)-1,1);
% for n=1:length(Nx)-1
%     legendInfo{n} = strcat('DNS\_',num2str(Nx(1)),'\times',num2str(Ny(1)),' vs EDFM\_',num2str(Nx(n+1)),'\times',num2str(Ny(n+1)));
%     plot( (1:t)*dt , Error_S(:,n) ,'LineWidth',2 );
%     hold on;
% end
% title('Saturation Error');
% legend(legendInfo);
% xlabel('Simulation Time [day], (\Deltat = 10^{-4} day)');
% ylabel('Error(S)');
% %% Manual Plot for 3 comparisons
% figure(2);
% hold on;
% plot( (1:t)*dt , Error_P(:,1) , '-'  , 'LineWidth',2 );
% plot( (1:t)*dt , Error_P(:,2) , '--' , 'LineWidth',2 );
% plot( (1:t)*dt , Error_P(:,3) , '-.' , 'LineWidth',2 );
% title('Pressure Error');
% legend('DNS(225\times225) vs EDFM(45\times45)' , 'DNS(225\times225) vs EDFM(25\times25)' , 'DNS(225\times225) vs EDFM(15\times15)');
% xlabel('Simulation Time [day], (\Deltat = 10^{-4} day)');
% ylabel('Error(P)');
% 
% %% Manual Plot for 3 comparisons
% figure(3);
% hold on;
% plot( (1:t)*dt , Error_S(:,1) , '-'  , 'LineWidth',2 );
% plot( (1:t)*dt , Error_S(:,2) , '--' , 'LineWidth',2 );
% plot( (1:t)*dt , Error_S(:,3) , '-.' , 'LineWidth',2 );
% title('Saturation Error');
% legend('DNS(225\times225) vs EDFM(45\times45)' , 'DNS(225\times225) vs EDFM(25\times25)' , 'DNS(225\times225) vs EDFM(15\times15)');
% xlabel('Simulation Time [day], (\Deltat = 10^{-4} day)');
% ylabel('Error(S)');
