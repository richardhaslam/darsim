% ERROR COMPARISON 
%DS = 'DS02';
%dsi = 1;
Problem = 'DARSIMRuns';
%% Read data from files
%Pressure
file = fopen(strcat('../Input/',Problem,'/FineScale/Output/FIMPressure.txt'), 'r');
X = fscanf(file, '%f %f %f %f %f %f %f %f %f %f', [10 Inf]);
FineGrid.P = X';
fclose(file);
file = fopen(strcat('../Input/',Problem,'/Constant/',DS,'/Output/Pressure.txt'), 'r');
X = fscanf(file, '%f %f %f %f %f %f %f %f %f %f', [10 Inf]);
Const2L.P = X';
fclose(file);
file = fopen(strcat('../Input/',Problem,'/Bilinear/', DS, '/Output/Pressure.txt'), 'r');
X = fscanf(file, '%f %f %f %f %f %f %f %f %f %f', [10 Inf]);
Bilin2L.P = X';
fclose(file);
% file = fopen(strcat('../Input/DARSim1/ADM/', DS,'/Constant_3Levels/Pressure.txt'), 'r');
% X = fscanf(file, '%f %f %f %f %f %f %f %f %f %f', [10 Inf]);
% Const3L.P = X';
% fclose(file);
% file = fopen(strcat('../Input/DARSim1/ADM/', DS, '/Bilinear_3Levels/Pressure.txt'), 'r');
% X = fscanf(file, '%f %f %f %f %f %f %f %f %f %f', [10 Inf]);
% Bilin3L.P = X';
% fclose(file);
%Saturation
file = fopen(strcat('../Input/',Problem,'/FineScale/Output/FIMSaturation.txt'), 'r');
X = fscanf(file, '%f %f %f %f %f %f %f %f %f %f', [10 Inf]);
FineGrid.S = X';
fclose(file);
file = fopen(strcat('../Input/',Problem,'/Constant/', DS, '/Output/Saturation.txt'), 'r');
X = fscanf(file, '%f %f %f %f %f %f %f %f %f %f', [10 Inf]);
Const2L.S = X';
fclose(file);
file = fopen(strcat('../Input/',Problem,'/Bilinear/', DS, '/Output/Saturation.txt'), 'r');
X = fscanf(file, '%f %f %f %f %f %f %f %f %f %f', [10 Inf]);
Bilin2L.S = X';
fclose(file);
% file = fopen(strcat('../Input/DARSim1/ADM/', DS, '/Constant_3Levels/Saturation.txt'), 'r');
% X = fscanf(file, '%f %f %f %f %f %f %f %f %f %f', [10 Inf]);
% Const3L.S = X';
% fclose(file);
% file = fopen(strcat('../Input/DARSim1/ADM/', DS, '/Bilinear_3Levels/Saturation.txt'), 'r');
% X = fscanf(file, '%f %f %f %f %f %f %f %f %f %f', [10 Inf]);
% Bilin3L.S = X';
% fclose(file);
%Stats
file = strcat('../Input/',Problem,'/FineScale/Output/FIMStat.txt');
FineGrid.Stats = load(file);
file = strcat('../Input/',Problem,'/Constant/', DS, '/Output/DLGRStat.txt');
Const2L.Stats = load(file);
file = strcat('../Input/',Problem,'/Bilinear/', DS, '/Output/DLGRStat.txt');
Bilin2L.Stats = load(file);
% file = strcat('../Input/DARSim1/ADM/', DS, '/Constant_3Levels/DLGRStat.txt');
% Const3L.Stats = load(file);
% file = strcat('../Input/DARSim1/ADM/', DS, '/Bilinear_3Levels/DLGRStat.txt');
% Bilin3L.Stats = load(file);

%Compute Error
for i = 1:10
    Const2L.errorP(i) =  norm(Const2L.P(:,i) -  FineGrid.P(:,i))./norm(FineGrid.P(:,i));
    Bilin2L.errorP(i) =  norm(Bilin2L.P(:,i) -  FineGrid.P(:,i))./norm(FineGrid.P(:,i));
%     Const3L.errorP(i) =  norm(Const3L.P(:,i) -  FineGrid.P(:,i))./norm(FineGrid.P(:,i));
%     Bilin3L.errorP(i) =  norm(Bilin3L.P(:,i) -  FineGrid.P(:,i))./norm(FineGrid.P(:,i));
    Const2L.errorS(i) =  norm(Const2L.S(:,i) -  FineGrid.S(:,i));%./norm(FineGrid.S(:,i), inf);
    Bilin2L.errorS(i) =  norm(Bilin2L.S(:,i) -  FineGrid.S(:,i));%./norm(FineGrid.S(:,i), inf);
%     Const3L.errorS(i) =  norm(Const3L.S(:,i) -  FineGrid.S(:,i));%./norm(FineGrid.S(:,i), inf);
%     Bilin3L.errorS(i) =  norm(Bilin3L.S(:,i) -  FineGrid.S(:,i));%./norm(FineGrid.S(:,i), inf);
end

%% Compute % of fine cells
FineCells = 11664;
Const2L.Stats(:, 4) = sum(Const2L.Stats(:, 4:6),2)./FineCells.*100; 
Bilin2L.Stats(:, 4) = sum(Bilin2L.Stats(:, 4:6),2)./FineCells.*100;
% Const3L.Stats(:, 4) = Const3L.Stats(:, 4)./FineCells.*100;
% Bilin3L.Stats(:, 4) = Bilin3L.Stats(:, 4)./FineCells.*100;

%% Plot % of active nodes
figure(150 + dsi)
axes('Position',[.05 .005 .90 .99],'xtick',[],'ytick',[],'box','on','handlevisibility','off')
hold on
title([strcat(DS,': ') ' Active nodes']);
plot(Const2L.Stats(:,1),Const2L.Stats(:,4), 'LineWidth', 3, 'color', 'red');
plot(Bilin2L.Stats(:,1),Bilin2L.Stats(:,4), 'LineWidth', 3, 'color', 'blue');
% plot(Const3L.Stats(:,1),Const3L.Stats(:,4), 'LineWidth', 3, 'color', 'green');
% plot(Bilin3L.Stats(:,1),Bilin3L.Stats(:,4), 'LineWidth', 3, 'color', 'yellow');
xlabel('Timestep');
ylabel('% of active nodes');
field1 = ['Const_{2L} - Av. = '  num2str(mean(Const2L.Stats(:,4)),4) ' %'];
field2 = ['Bilin_{2L} - Av. = '  num2str(mean(Bilin2L.Stats(:,4)),4) ' %'];
% field3 = ['Const_{3L} - Av. = '  num2str(mean(Const3L.Stats(:,4)),4) ' %'];
% field4 = ['Bilin_{3L} - Av. = '  num2str(mean(Bilin3L.Stats(:,4)),4) ' %'];
legend({field1, field2}, 'FontSize', 32, 'FontWeight', 'bold', 'Location', 'best');
set(gca,'fontsize',28);

%% Print Data file for pgfplots
file = fopen(strcat('../Input/', Problem,'/', Problem,'Const2L_AN_', DS, '.dat'), 'w');
fprintf(file, '%10.5f %10.5f\n', [Const2L.Stats(:,1), Const2L.Stats(:,4)]');
fclose(file);
file = fopen(strcat('../Input/', Problem, '/', Problem,'Bilin2L_AN_', DS, '.dat'), 'w');
fprintf(file, '%10.5f %10.5f\n', [Bilin2L.Stats(:,1), Bilin2L.Stats(:,4)]');
fclose(file);
% file = fopen(strcat('../Input/DARSim1/ADM/', DS, '/Constant_3Levels/DARSim1Const3L_AN_', DS, '.dat'), 'w');
% fprintf(file, '%10.5f %10.5f\n', [Const3L.Stats(:,1), Const3L.Stats(:,4)]');
% fclose(file);
% file = fopen(strcat('../Input/DARSim1/ADM/', DS, '/Bilinear_3Levels/DARSim1Bilin3L_AN_', DS, '.dat'), 'w');
% fprintf(file, '%10.5f %10.5f\n', [Bilin3L.Stats(:,1), Bilin3L.Stats(:,4)]');
% fclose(file);

%% Number of nonlinear iterations
Ref = sum(FineGrid.Stats(:,3))
Const2L_nl = sum(Const2L.Stats(:,3))
Bilin2L_nl = sum(Bilin2L.Stats(:,3))
% Const3L_nl = sum(Const3L.Stats(:,3));
% Bilin3L_nl = sum(Bilin3L.Stats(:,3));

%% Error
Const2L_EP = mean(Const2L.errorP);
Bilin2L_EP = mean(Bilin2L.errorP);
% Const3L_EP = mean(Const3L.errorP)*100;
% Bilin3L_EP = mean(Bilin3L.errorP)*100;
Const2L_ES = mean(Const2L.errorS)/FineCells;
Bilin2L_ES = mean(Bilin2L.errorS)/FineCells;
% Const3L_ES = mean(Const3L.errorS);
% Bilin3L_ES = mean(Bilin3L.errorS);

%% Write them in a file
file = fopen(strcat('../Input/',Problem, '/',Problem,'_', DS, '.dat'), 'w');
fprintf(file, '%20s %10s %10s \n','', 'Const2L', 'Bilin2L');%, 'Const3L', 'Bilin3L');
fprintf(file, '%20s %10.2f %10.2f \n', '% active nodes', [mean(Const2L.Stats(:,4)), mean(Bilin2L.Stats(:,4))]);%, mean(Const3L.Stats(:,4)), mean(Bilin3L.Stats(:,4))]);
fprintf(file, '%20s %10.0f %10.0f \n', 'NnLinear', [Const2L_nl, Bilin2L_nl]);%, Const3L_nl, Bilin3L_nl]);
fprintf(file, '%20s %10.5f %10.5f \n', 'Pressure error', [Const2L_EP, Bilin2L_EP]);%, Const3L_EP, Bilin3L_EP]);
fprintf(file, '%20s %10.5f %10.5f \n', 'Saturation error', [Const2L_ES, Bilin2L_ES]);%, Const3L_ES, Bilin3L_ES]);
fclose(file);