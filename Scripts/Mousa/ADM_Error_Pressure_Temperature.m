close all; clear all; clc;
Directory = 'D:\SURFdrive\Simulation\Results_Comparison\Results_Comparison_pEDFM_ADM\SinglePhaseGeo_2D_136x136_30frac_pEDFM\MixedConductiveFractures\Heterogeneous\Input_files';

Nx = 136;
Ny = 136;
Nz = 001;
Nf = 1044;
N_tstep = 11;
N_Comparison = 4;
tol = [5,10,20,50];

Error_P = zeros(N_tstep,N_Comparison);
Error_T = zeros(N_tstep,N_Comparison);
ADMStats = cell(N_Comparison,1);

%% Looping over time to compute pressure and saturation errors
for t = 1:N_tstep
    
    fprintf('Checking timestep: %3.0f\n', t);
	
    % Reading finescale files
	File = strcat(Directory, '\FineScale\SinglePhaseGeo_Sol',num2str(t),'.txt');
	FileID = fopen(File, 'r');
	Matrix = textscan(FileID, '%s');
	Matrix = str2double(Matrix{1,1});
	Matrix = reshape(Matrix,[5,Nx*Ny*Nz+Nf])';
	fclose(FileID);
	P_FineScale = Matrix(:,2);
	T_FineScale = Matrix(:,4);
	
	% Reading ADM files
    for n = 1:N_Comparison
        File = strcat(Directory, '\ADM_',num2str(tol(n)),'\SinglePhaseGeo_Sol',num2str(t),'.txt');
		FileID = fopen(File, 'r');
		Matrix = textscan(FileID, '%s');
		Matrix = str2double(Matrix{1,1});
		Matrix = reshape(Matrix,[5,Nx*Ny*Nz+Nf])';
		fclose(FileID);
		P_ADM{n} = Matrix(:,2);
		T_ADM{n} = Matrix(:,4);
        Error_P(t,n) = norm( (P_FineScale - P_ADM{n}) ) / norm( P_FineScale );
        Error_T(t,n) = norm( (T_FineScale - T_ADM{n}) ) / norm( T_FineScale );
    end
end

%% Reading ADMStats files
for n = 1:N_Comparison
    File = strcat(Directory, '\ADM_',num2str(tol(n)),'\ADMStats.txt');
    FileID = fopen(File, 'r');
    Matrix = textscan(FileID, '%s');
    Matrix = str2double(Matrix{1,1});
    Matrix = reshape(Matrix,[4,length(Matrix)/4])';
    fclose(FileID);
    ADMStats{n} = Matrix(:,4);
end

%% Computing average values
Error_P_Ave = zeros(1,length(tol));
Error_T_Ave = zeros(1,length(tol));
ADMStats_Ave = zeros(1,length(tol));
for n = 1 : N_Comparison
    Error_P_Ave(n) = mean(Error_P(:,n));
    Error_T_Ave(n) = mean(Error_T(:,n));
    ADMStats_Ave(n) = mean(ADMStats{n});
end

%% Saving the Workspace
save(strcat(Directory,'\matlab.mat'));

%% Plotting attributes
figLineType = { '-' , '--' , '-.' , '-*' };
legendText = cell(N_Comparison,1);
for n = 1:N_Comparison
    legendText(n) = strcat( {'tol '} , char(916) , 'T=' , num2str(tol(n)) );
end
%% Plotting pressure error
fig = figure;
hold on;
for n = 1:N_Comparison
    plot( (1:t) , Error_P(:,n) , figLineType{n}  , 'LineWidth',2 );
end
xlim([1,N_tstep]);
set(gca,'fontsize',15);
set(gcf, 'Position', [100, 100, 800, 500]);
grid on;
title('Pressure Error');
legend(legendText,'FontSize', 15, 'Location', 'north');
xlabel('Time-step','FontSize',15);
ylabel('Error(P)','FontSize',15);
saveas(fig, strcat(Directory, '\Pressure_Error.fig'));
saveas(fig, strcat(Directory, '\Pressure_Error.png'));

%% Plotting temperature error
fig = figure;
hold on;
for n = 1:N_Comparison
    plot( (1:t) , Error_T(:,n) , figLineType{n}  , 'LineWidth',2 );
end
xlim([1,N_tstep]);
set(gca,'fontsize',15);
set(gcf, 'Position', [100, 100, 800, 500]);
grid on;
title('Temperature Error');
legend(legendText,'FontSize', 15, 'Location', 'north');
xlabel('Time-step','FontSize',15);
ylabel('Error(S)','FontSize',15);
saveas(fig, strcat(Directory, '\Temperature_Error.fig'));
saveas(fig, strcat(Directory, '\Temperature_Error.png'));

%% Plotting active grid cells
fig = figure;
hold on;
for n = 1:N_Comparison
    plot( (1:N_tstep) , ADMStats{n}(1:N_tstep) , figLineType{n}  , 'LineWidth',2 );
end
xlim([1,N_tstep]);
set(gca,'fontsize',15);
set(gcf, 'Position', [100, 100, 800, 500]);
grid on;
title('Percentage of Active Grid Cells over Time-steps');
legend(legendText,'FontSize', 15, 'Location', 'north');
xlabel('Time-step','FontSize',15);
ylabel('% Active Grid Cells','FontSize',15);
saveas(fig, strcat(Directory, '\Active_Grids.fig'));
saveas(fig, strcat(Directory, '\Active_Grids.png'));

%% Plotting average values for pressure error, temperature error and active grid cells for different tolerances
fig = figure;
hold on;
yyaxis left; plot( tol , Error_P_Ave , '--ob'  , 'LineWidth',2 , 'MarkerSize',12 );
yyaxis right; plot( tol , Error_T_Ave , '--*r' , 'LineWidth',2 , 'MarkerSize',12 );
yyaxis right; plot( tol , ADMStats_Ave/100 , '--xg' , 'LineWidth',2 , 'MarkerSize',12 );
set(gca,'fontsize',15);
set(gcf, 'Position', [100, 100, 800, 500]);
grid on;
legend({'Average Perssure Error (left axis)', 'Average Temperature Error (right axis)', 'Avergare Active Grid Cells [-] (right axis)'},'FontSize', 15, 'Location', 'north');
xlabel('ADM Tolerance','FontSize',15);
ylabel('[-]','FontSize',15);
saveas(fig, strcat(Directory, '\Average_Values.fig'));
saveas(fig, strcat(Directory, '\Average_Values.png'));