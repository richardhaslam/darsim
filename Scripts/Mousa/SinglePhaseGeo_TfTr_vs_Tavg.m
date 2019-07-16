close all; clear all; clc;
Directory = 'D:\SURFdrive\Simulation\Results_Comparison\Results_Comparison_FGADM\TestCase2_FineScale_TfTr_vs_Tavg\Input\';

Nx = 81;
Ny = 81;
Nz = 001;
Nf = 1188;%1665;
N_tstep = 25;
T_Inj = 300;
T_Res = 400;

Error_P  = zeros(N_tstep,1);
Error_Tf = zeros(N_tstep,1);
Error_Tr = zeros(N_tstep,1);
TfTr_Diff = zeros(N_tstep,1);

%% Looping over time
for t = 1:N_tstep
    
    fprintf('Checking timestep: %3.0f\n', t);
	
    % Reading Finescale Files
	File = strcat(Directory, 'SinglePhaseGeo_FineScale_TfTr_Sol',num2str(t),'.txt');
	FileID = fopen(File, 'r');
	Matrix = textscan(FileID, '%s');
	Matrix = str2double(Matrix{1,1});
	Matrix = reshape(Matrix,[5,Nx*Ny*Nz+Nf])';
	fclose(FileID);
	P_TfTr  = Matrix(:,2);
	Tf_TfTr = Matrix(:,4);
    Tr_TfTr = Matrix(1:Nx*Ny*Nz,5);
	
    File = strcat(Directory, 'SinglePhaseGeo_FineScale_Tavg_Sol', num2str(t),'.txt');
    FileID = fopen(File, 'r');
	Matrix = textscan(FileID, '%s');
	Matrix = str2double(Matrix{1,1});
	Matrix = reshape(Matrix,[5,Nx*Ny*Nz+Nf])';
	fclose(FileID);
	P_Tavg  = Matrix(:,2);
	Tf_Tavg = Matrix(:,4);
    Tr_Tavg = Matrix(1:Nx*Ny*Nz,5);
    Error_P (t,1) = norm( P_TfTr  - P_Tavg  ) / norm( P_TfTr );
    Error_Tf(t,1) = norm( Tf_TfTr - Tf_Tavg ) / norm( Tf_TfTr );
    Error_Tr(t,1) = norm( Tr_TfTr - Tr_Tavg ) / norm( Tr_TfTr );
    TfTr_Diff(t,1) = max(abs(Tf_TfTr(1:Nx*Ny*Nz) - Tr_TfTr)) / abs(T_Inj - T_Res);
end


%% Plot for Pressure Error
figure;
hold on;
yyaxis right; plot( (1:t) , TfTr_Diff, '-r'  , 'LineWidth',2 );
yyaxis left;  plot( (1:t) , Error_P  , '--b' , 'LineWidth',2 );
yyaxis left;  plot( (1:t) , Error_Tf , '-.g' , 'LineWidth',2 );
xlim([1,N_tstep]);
set(gca,'fontsize',15);
set(gcf, 'Position', [100, 100, 800, 500]);
grid on;
legend({'$||p_{full} - p_{reduced}|| / ||p_{full}||$ (left axis)', '$||T_{f,full}- T_{avg,reduced}|| / ||T_{f,full}||$ (left axis)','$max(|T_{f} - T_{r}|) / (T_{max} - T_{min}) $ (right axis)'},'Interpreter','latex','FontSize', 15);
xlabel('Time-step','FontSize',15);

% %% Plot for Temperature
% figure;
% hold on;
% plot( (1:t) , Error_Tf(:,1) , '-'  , 'LineWidth',2 );
% plot( (1:t) , Error_Tf(:,2) , '--' , 'LineWidth',2 );
% plot( (1:t) , Error_Tf(:,3) , '-.' , 'LineWidth',2 );
% plot( (1:t) , Error_Tf(:,4) , '-*' , 'LineWidth',2 );
% xlim([1,N_tstep]);
% set(gca,'fontsize',15);
% set(gcf, 'Position', [100, 100, 800, 500]);
% grid on;
% title('Temperature Error');
% legend({'Tol 5', 'Tol 10', 'Tol 20', 'Tol 50'},'FontSize', 15);
% xlabel('Time-step','FontSize',15);
% ylabel('Error(Tf)','FontSize',15);

%% %% Plot for Rock Temperature
% figure;
% hold on;
% plot( (1:t) , Error_Tr(:,1) , '-'  , 'LineWidth',2 );
% plot( (1:t) , Error_Tr(:,2) , '--' , 'LineWidth',2 );
% plot( (1:t) , Error_Tr(:,3) , '-.' , 'LineWidth',2 );
% plot( (1:t) , Error_Tr(:,4) , '-*' , 'LineWidth',2 );
% xlim([1,N_tstep]);
% set(gca,'fontsize',15);
% set(gcf, 'Position', [100, 100, 800, 500]);
% grid on;
% title('Rock Temperature Error');
% legend({'Tol 0.1', 'Tol 0.3', 'Tol 0.5', 'Tol 0.8'},'FontSize', 15);
% xlabel('Time-step','FontSize',15);
% ylabel('Error(S)','FontSize',15);

fclose('all');