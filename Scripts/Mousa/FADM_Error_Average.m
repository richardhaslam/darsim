Directory = 'D:\SURFdrive\Simulation\Results_Comparison\Results_Comparison_FADM\F_ADM_2D_Coupled_3wells_hetro_patchy\Input_2D\';
%load Result_Comparison_FADM_2D_Error_Workspace_Coupled.mat;
N_tstep = 50; % define yourself
N_Comparison = 4; % define yourself
tol = [0.1,0.3,0.5,0.8];
ADMStats = cell(N_Comparison,1);
%% Reading ADMStats Files
for n = 1:N_Comparison
    File = strcat(Directory, 'ImmHomo_FADM_', num2str(tol(n)),'_ADMStats.txt');
    FileID = fopen(File, 'r');
    Matrix = textscan(FileID, '%s');
    Matrix = str2double(Matrix{1,1});
    Matrix = reshape(Matrix,[4,length(Matrix)/4])';
    fclose(FileID);
    ADMStats{n} = Matrix(:,4);
end

%% Plot for Active Grid Cells
figure;
hold on;
plot( (1:N_tstep) , ADMStats{1}(1:N_tstep) , '-'  , 'LineWidth',2 );
plot( (1:N_tstep) , ADMStats{2}(1:N_tstep) , '--' , 'LineWidth',2 );
plot( (1:N_tstep) , ADMStats{3}(1:N_tstep) , '-.' , 'LineWidth',2 );
plot( (1:N_tstep) , ADMStats{4}(1:N_tstep) , '-*' , 'LineWidth',2 );
xlim([1,N_tstep]);
set(gca,'fontsize',15);
set(gcf, 'Position', [100, 100, 800, 500]);
grid on;
%title('Percentage of Active Grid Cells over Time-steps');
legend({'Tol 0.1', 'Tol 0.3', 'Tol 0.5', 'Tol 0.8'},'FontSize', 15);
xlabel('Time-step','FontSize',15);
ylabel('% Active Grid Cells','FontSize',15);

%% Plot of Average Pressure Error, Saturation Error and Active Grid Cells versus Different Tolerances
Error_P_Ave = zeros(1,length(tol));
Error_S_Ave = zeros(1,length(tol));
ADMStats_Ave = zeros(1,length(tol));
for n = 1 : length(tol)
    Error_P_Ave(n) = mean(Error_P(:,n));
    Error_S_Ave(n) = mean(Error_S(:,n));
    ADMStats_Ave(n) = mean(ADMStats{n});
end

figure;
hold on;
yyaxis left; plot( tol , Error_P_Ave , '--or'  , 'LineWidth',2 , 'MarkerSize',12 );
yyaxis right; plot( tol , Error_S_Ave , '--*b' , 'LineWidth',2 , 'MarkerSize',12 );
yyaxis right; plot( tol , ADMStats_Ave/100 , '--xg' , 'LineWidth',2 , 'MarkerSize',12 );
set(gca,'fontsize',15);
set(gcf, 'Position', [100, 100, 800, 500]);
grid on;
% title('Average of Pressure Error, Saturation Error and Percentage of Active Grid Cells over the Whole Simulation Time');
legend({'Average Perssure Error (left axis)', 'Average Saturation Error (right axis)', 'Avergare Active Grid Cells [-] (right axis)'},'FontSize', 15);
xlabel('F-ADM Tolerance','FontSize',15);
ylabel('[-]','FontSize',15);