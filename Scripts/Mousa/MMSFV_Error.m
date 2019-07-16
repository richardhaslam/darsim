close all; clear all; clc;
Directory = 'D:/SURFdrive/Simulation/Results_Comparison/Results_Comparison_MMSFV/Input/';

N_TestCase = 6;
TestCase_Name = 'SinglePhase';

Lx = 75;  Nx = 135;  dx = Lx/Nx;
Ly = 75;  Ny = 135;  dy = Ly/Ny;
Lz = 01;  Nz = 001;  dz = Lz/Nz;
 
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
    end
    for n = 1:N_TestCase - 1
        Error_P(n,t) = norm( (P(:,1) - P(:,n+1)) ) / length(P(:,1));
    end
end