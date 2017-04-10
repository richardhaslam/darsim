% ERROR COMPARISON for Heterogeneous Problems
DS = 'DS02';
dsi = 1;
Problem = 'SPE10BRuns';
%% Read data from files
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
file = fopen(strcat('../Input/',Problem,'/MS/', DS, '/Output/Pressure.txt'), 'r');
X = fscanf(file, '%f %f %f %f %f %f %f %f %f %f', [10 Inf]);
MS2L.P = X';
fclose(file);
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
file = fopen(strcat('../Input/',Problem,'/MS/', DS, '/Output/Saturation.txt'), 'r');
X = fscanf(file, '%f %f %f %f %f %f %f %f %f %f', [10 Inf]);
MS2L.S = X';
fclose(file);
%Stats
file = strcat('../Input/',Problem,'/FineScale/Output/FIMStat.txt');
FineGrid.Stats = load(file);
file = strcat('../Input/',Problem,'/Constant/', DS, '/Output/DLGRStat.txt');
Const2L.Stats = load(file);
file = strcat('../Input/',Problem,'/Bilinear/', DS, '/Output/DLGRStat.txt');
Bilin2L.Stats = load(file);
file = strcat('../Input/',Problem,'/MS/', DS, '/Output/DLGRStat.txt');
MS2L.Stats = load(file);


%Compute Error
for i = 1:10
    Const2L.errorP(i) =  norm(Const2L.P(:,i) -  FineGrid.P(:,i))./norm(FineGrid.P(:,i));
    Bilin2L.errorP(i) =  norm(Bilin2L.P(:,i) -  FineGrid.P(:,i))./norm(FineGrid.P(:,i));
    MS2L.errorP(i) = norm(MS2L.P(:,i) -  FineGrid.P(:,i))./norm(FineGrid.P(:,i));
    Const2L.errorS(i) =  norm(Const2L.S(:,i) -  FineGrid.S(:,i));
    Bilin2L.errorS(i) =  norm(Bilin2L.S(:,i) -  FineGrid.S(:,i));
    MS2L.errorS(i) =  norm(MS2L.S(:,i) -  FineGrid.S(:,i));
end

%% Compute % of fine cells
%% Compute % of fine cells
FineCells = 11664;
Const2L.Stats(:, 4) = sum(Const2L.Stats(:, 4:6),2)./FineCells.*100; 
Bilin2L.Stats(:, 4) = sum(Bilin2L.Stats(:, 4:6),2)./FineCells.*100;
MS2L.Stats(:, 4) = sum(MS2L.Stats(:, 4:6),2)./FineCells.*100;
%% Plot % of active nodes
figure(150 + dsi)
axes('Position',[.05 .005 .90 .99],'xtick',[],'ytick',[],'box','on','handlevisibility','off')
hold on
title([strcat(DS,': ') ' Active nodes']);
plot(Const2L.Stats(:,1),Const2L.Stats(:,4), 'LineWidth', 3, 'color', 'red');
plot(Bilin2L.Stats(:,1),Bilin2L.Stats(:,4), 'LineWidth', 3, 'color', 'blue');
xlabel('Timestep');
ylabel('% of active nodes');
field1 = ['Const_{2L} - Av. = '  num2str(mean(Const2L.Stats(:,4)),4) ' %'];
field2 = ['Bilin_{2L} - Av. = '  num2str(mean(Bilin2L.Stats(:,4)),4) ' %'];
legend({field1, field2}, 'FontSize', 32, 'FontWeight', 'bold', 'Location', 'best');
set(gca,'fontsize',28);

%% Active nodes: print Data file for pgfplots
%Constant 2L
file = fopen(strcat('../Input/', Problem,'/', Problem,'Const2L_AN_', DS, '.dat'), 'w');
fprintf(file, '%10.5f %10.5f\n', [Const2L.Stats(:,1), Const2L.Stats(:,4)]');
fclose(file);
file = fopen(strcat('../Input/', Problem, '/', Problem,'Bilin2L_AN_', DS, '.dat'), 'w');
fprintf(file, '%10.5f %10.5f\n', [Bilin2L.Stats(:,1), Bilin2L.Stats(:,4)]');
fclose(file);
file = fopen(strcat('../Input/', Problem, '/',Problem,'MS2L_AN_', DS, '.dat'), 'w');
fprintf(file, '%10.5f %10.5f\n', [MS2L.Stats(:,1), MS2L.Stats(:,4)]');
fclose(file);

%% Number of nonlinear iterations
Ref = sum(FineGrid.Stats(:,3));
Const2L_nl = sum(Const2L.Stats(:,3));
Bilin2L_nl = sum(Bilin2L.Stats(:,3));
MS2L_nl = sum(MS2L.Stats(:,3));

%% Error
Const2L_EP = mean(Const2L.errorP);
Bilin2L_EP = mean(Bilin2L.errorP);
MS2L_EP = mean(MS2L.errorP);
figure
pcolor(reshape(abs(Bilin2L.P(:,5) - FineGrid.P(:,5)), 216, 54)');
colormap(jet);
colorbar;
axis('image');
Const2L_ES = mean(Const2L.errorS)/FineCells;
Bilin2L_ES = mean(Bilin2L.errorS)/FineCells;
MS2L_ES = mean(MS2L.errorS)/FineCells;

%% Write them in a file
file = fopen(strcat('../Input/',Problem, '/',Problem,'_', DS, '.dat'), 'w');
fprintf(file, '%20s %10s %10s %10s \n','', 'Const2L', 'Bilin2L', 'MS2L');
fprintf(file, '%20s %10.2f %10.2f %10.2f \n', '% active nodes', [mean(Const2L.Stats(:,4)), mean(Bilin2L.Stats(:,4)), mean(MS2L.Stats(:,4))]);
fprintf(file, '%20s %10.0f %10.0f %10.0f \n', 'NnLinear', [Const2L_nl, Bilin2L_nl, MS2L_nl]);
fprintf(file, '%20s %10.5f %10.5f %10.5f \n', 'Pressure error', [Const2L_EP, Bilin2L_EP, MS2L_EP]);
fprintf(file, '%20s %10.5f %10.5f %10.5f \n', 'Saturation error', [Const2L_ES, Bilin2L_ES, MS2L_ES]);
fclose(file);