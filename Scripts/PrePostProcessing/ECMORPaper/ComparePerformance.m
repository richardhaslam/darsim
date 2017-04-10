%PERFORMANCE COMPARISON

%% Read stats from files
file = '../ResultsAndImages/FIM/Homogeneous/FIMStat.txt';
FineGrid = load(file);
file = '../ResultsAndImages/ADM/DS003/Homogeneous/Constant_2Levels/DLGRStat.txt';
Const2L = load(file);
file = '../ResultsAndImages/ADM/DS003/Homogeneous/Bilinear_2Levels/DLGRStat.txt';
Bilin2L = load(file);
file = '../ResultsAndImages/ADM/DS003/Homogeneous/Constant_3Levels/DLGRStat.txt';
Const3L = load(file);
file = '../ResultsAndImages/ADM/DS003/Homogeneous/Bilinear_3Levels/DLGRStat.txt';
Bilin3L = load(file);

%% Compute % of fine cells
FineCells = 11664;
Const2L(:, 4) =  Const2L(:, 4)./FineCells.*100; 
Bilin2L(:, 4) = Bilin2L(:, 4)./FineCells.*100;
Const3L(:, 4) = Const3L(:, 4)./FineCells.*100;
Bilin3L(:, 4) = Bilin3L(:, 4)./FineCells.*100;

%% Plot % of active nodes
figure(150)
axes('Position',[.05 .005 .90 .99],'xtick',[],'ytick',[],'box','on','handlevisibility','off')
hold on
title('Active nodes');
plot(Const2L(:,1),Const2L(:,4), 'LineWidth', 3, 'color', 'red');
plot(Bilin2L(:,1),Bilin2L(:,4), 'LineWidth', 3, 'color', 'blue');
plot(Const3L(:,1),Const3L(:,4), 'LineWidth', 3, 'color', 'green');
plot(Bilin3L(:,1),Bilin3L(:,4), 'LineWidth', 3, 'color', 'yellow');
xlabel('Timestep');
ylabel('% of active nodes');
field1 = ['Const_{2L} - Av. = '  num2str(mean(Const2L(:,4)),4) ' %'];
field2 = ['Bilin_{2L} - Av. = '  num2str(mean(Bilin2L(:,4)),4) ' %'];
field3 = ['Const_{3L} - Av. = '  num2str(mean(Const3L(:,4)),4) ' %'];
field4 = ['Bilin_{3L} - Av. = '  num2str(mean(Bilin3L(:,4)),4) ' %'];
legend({field1, field2, field3, field4}, 'FontSize', 32, 'FontWeight', 'bold', 'Location', 'best');
set(gca,'fontsize',28);
%% Number of nonlinear iterations
Ref = sum(FineGrid(:,3))
Const2L_nl = sum(Const2L(:,3))
Bilin2L_nl = sum(Bilin2L(:,3))
Const3L_nl = sum(Const3L(:,3))
Bilin3L_nl = sum(Bilin3L(:,3))

%% Read stats from files
file = '../ResultsAndImages/FIM/Homogeneous/FIMTimings.txt';
FineGrid = load(file);
file = '../ResultsAndImages/ADM/DS003/Homogeneous/Constant_2Levels/DLGRTimings.txt';
Const2L = load(file);
file = '../ResultsAndImages/ADM/DS003/Homogeneous/Bilinear_2Levels/DLGRTimings.txt';
Bilin2L = load(file);
file = '../ResultsAndImages/ADM/DS003/Homogeneous/Constant_3Levels/DLGRTimings.txt';
Const3L = load(file);
file = '../ResultsAndImages/ADM/DS003/Homogeneous/Bilinear_3Levels/DLGRTimings.txt';
Bilin3L = load(file);
