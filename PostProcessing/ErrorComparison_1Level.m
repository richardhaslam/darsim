% ERROR for Heterogeneous problems for 1Level ADM
%% Read data from files
file = fopen('../ResultsAndImages/FIM/SPE10T/FIMPressure.txt', 'r');
X = fscanf(file, '%f %f %f %f %f %f %f %f %f %f', [10 Inf]);
FineGrid.P = X';
fclose(file);
file = fopen('../ResultsAndImages/FIM/SPE10T/FIMSaturation.txt', 'r');
X = fscanf(file, '%f %f %f %f %f %f %f %f %f %f', [10 Inf]);
FineGrid.S = X';
fclose(file);
file = fopen('../ResultsAndImages/ADM/DS01/SPE10T/MS_1Level/Pressure.txt', 'r');
X = fscanf(file, '%f %f %f %f %f %f %f %f %f %f', [10 Inf]);
MS1L.P = X';
fclose(file);
file = fopen('../ResultsAndImages/ADM/DS01/SPE10T/MS_1Level/Saturation.txt', 'r');
X = fscanf(file, '%f %f %f %f %f %f %f %f %f %f', [10 Inf]);
MS1L.S = X';
fclose(file);
file = fopen('../ResultsAndImages/ADM/DS01/SPE10T/Constant_1Level/Pressure.txt', 'r');
X = fscanf(file, '%f %f %f %f %f %f %f %f %f %f', [10 Inf]);
Const1L.P = X';
fclose(file);
file = fopen('../ResultsAndImages/ADM/DS01/SPE10T/Constant_1Level/Saturation.txt', 'r');
X = fscanf(file, '%f %f %f %f %f %f %f %f %f %f', [10 Inf]);
Const1L.S = X';
fclose(file);
%Stats
file = '../ResultsAndImages/FIM/SPE10T/FIMStat.txt';
FineGrid.Stats = load(file);
file = '../ResultsAndImages/ADM/DS01/SPE10T/Constant_1Level/DLGRStat.txt';
Const1L.Stats = load(file);
file = '../ResultsAndImages/ADM/DS01/SPE10T/MS_1Level/DLGRStat.txt';
MS1L.Stats = load(file);

%% Compute Error
for i = 1:10
    Const1L.errorP(i) =  norm(Const1L.P(:,i) -  FineGrid.P(:,i))./norm(FineGrid.P(:,i));
    MS1L.errorP(i) =  norm(MS1L.P(:,i) -  FineGrid.P(:,i))./norm(FineGrid.P(:,i));
    Const1L.errorS(i) = norm(MS1L.S(:,i) -  FineGrid.S(:,i));%./norm(FineGrid.S(:,i));
    MS1L.errorS(i) =  norm(MS1L.S(:,i) -  FineGrid.S(:,i));%./norm(FineGrid.S(:,i));
end

%% Compute % of fine cells
FineCells = 11664;
Const1L.Stats(:, 4) = Const1L.Stats(:, 4)./FineCells.*100; 
MS1L.Stats(:, 4) = MS1L.Stats(:, 4)./FineCells.*100;

figure(150)
axes('Position',[.05 .005 .90 .99],'xtick',[],'ytick',[],'box','on','handlevisibility','off')
hold on
title('Active nodes');
plot(Const1L.Stats(:,1),Const1L.Stats(:,4), 'LineWidth', 3, 'color', 'red');
plot(MS1L.Stats(:,1),MS1L.Stats(:,4), 'LineWidth', 3, 'color', 'blue');
xlabel('Timestep');
ylabel('% of active nodes');
field1 = ['Const_{1L} - Av. = '  num2str(mean(Const1L.Stats(:,4)),4) ' %'];
field2 = ['MS_{1L} - Av. = '  num2str(mean(MS1L.Stats(:,4)),4) ' %'];
legend({field1, field2}, 'FontSize', 32, 'FontWeight', 'bold', 'Location', 'best');
set(gca,'fontsize',28);

%% Number of nonlinear iterations
disp('number of nonlinear iterations');
Ref_nl = sum(FineGrid.Stats(:,3))
Const1L_nl = sum(Const1L.Stats(:,3))
MS1L_nl = sum(MS1L.Stats(:,3))

%% Error
disp('Pressure error');
Const1L_EP = mean(Const1L.errorP)*100
MS1L_EP = mean(MS1L.errorP)*100
disp('Saturation error')
Const1L_ES = mean(Const1L.errorS)
MS1L_ES = mean(MS1L.errorS)