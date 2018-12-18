close all; clear all; clc;
Directory = 'Input_2D/';

Nx = 81;
Ny = 81;
Nz = 1;
Nf = 0;
N_tstep = 50;
N_Comparison = 4;
tol = [0.1,0.3,0.5,0.8];

Error_P = zeros(N_tstep,N_Comparison);
Error_S = zeros(N_tstep,N_Comparison);

%% Looping over time
for t = 1:N_tstep
    
    fprintf('Checking timestep: %3.0f\n', t);
	
    % Reading Finescale Files
	File = strcat(Directory, 'ImmHomo_FineScale_Sol',num2str(t),'.txt');
	FileID = fopen(File, 'r');
	Matrix = textscan(FileID, '%s');
	Matrix = str2double(Matrix{1,1});
	Matrix = reshape(Matrix,[3,Nx*Ny*Nz+Nf])';
	fclose(FileID);
	P_FineScale = Matrix(:,2);
	S_FineScale = Matrix(:,3);
	
	% Reading F-ADM Files
    for n = 1:N_Comparison
        File = strcat(Directory, 'ImmHomo_ADM_', num2str(tol(n)),'_Sol',num2str(t),'.txt');
		FileID = fopen(File, 'r');
		Matrix = textscan(FileID, '%s');
		Matrix = str2double(Matrix{1,1});
		Matrix = reshape(Matrix,[3,Nx*Ny*Nz+Nf])';
		fclose(FileID);
		P_FADM{n} = Matrix(:,2);
		S_FADM{n} = Matrix(:,3);
        Error_P(t,n) = norm( (P_FineScale - P_FADM{n}) ) / norm( P_FineScale );
        Error_S(t,n) = norm( (S_FineScale - S_FADM{n}) ) / norm( S_FineScale );
    end
end

%% Plot for Pressure Error
figure;
hold on;
plot( (1:t) , Error_P(:,1) , '-'  , 'LineWidth',2 );
plot( (1:t) , Error_P(:,2) , '--' , 'LineWidth',2 );
plot( (1:t) , Error_P(:,3) , '-.' , 'LineWidth',2 );
plot( (1:t) , Error_P(:,4) , '-*' , 'LineWidth',2 );
xlim([1,N_tstep]);
set(gca,'fontsize',15);
set(gcf, 'Position', [100, 100, 800, 500]);
grid on;
title('Pressure Error');
legend({'Tol 0.1', 'Tol 0.3', 'Tol 0.5', 'Tol 0.8'},'FontSize', 15);
xlabel('Time-step','FontSize',15);
ylabel('Error(P)','FontSize',15);

%% Plot for Saturation Error
figure;
hold on;
plot( (1:t) , Error_S(:,1) , '-'  , 'LineWidth',2 );
plot( (1:t) , Error_S(:,2) , '--' , 'LineWidth',2 );
plot( (1:t) , Error_S(:,3) , '-.' , 'LineWidth',2 );
plot( (1:t) , Error_S(:,4) , '-*' , 'LineWidth',2 );
xlim([1,N_tstep]);
set(gca,'fontsize',15);
set(gcf, 'Position', [100, 100, 800, 500]);
grid on;
title('Saturation Error');
legend({'Tol 0.1', 'Tol 0.3', 'Tol 0.5', 'Tol 0.8'},'FontSize', 15);
xlabel('Time-step','FontSize',15);
ylabel('Error(S)','FontSize',15);