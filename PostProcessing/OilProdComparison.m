%%%%%% Oil Production Comparison %%%%%%%
clear all
close all

%% Read stats from files
file = '../ResultsAndImages/FIM/Homogeneous/FIMOilProd.txt';
FineGrid = load(file);
file = '../ResultsAndImages/ADM/DS003/Homogeneous/Constant_2Levels/DLGROilProd.txt';
Const2L = load(file);
file = '../ResultsAndImages/ADM/DS003/Homogeneous/Bilinear_2Levels/DLGROilProd.txt';
Bilin2L = load(file);
file = '../ResultsAndImages/ADM/DS003/Homogeneous/Constant_3Levels/DLGROilProd.txt';
Const3L = load(file);
file = '../ResultsAndImages/ADM/DS003/Homogeneous/Bilinear_3Levels/DLGROilProd.txt';
Bilin3L = load(file);

%% Plot
figure(500)
axes('Position',[.05 .005 .90 .99],'xtick',[],'ytick',[],'box','on','handlevisibility','off');
hold on
plot(FineGrid(:,1), FineGrid(:,2),'LineWidth', 3, 'color', 'black');
plot(Const2L(:,1), Const2L(:,2),'LineWidth', 3, 'color', 'red');
plot(Bilin2L(:,1), Bilin2L(:,2),'LineWidth', 3, 'color', 'blue');
plot(Const3L(:,1), Const3L(:,2),'LineWidth', 3, 'color', 'green');
plot(Bilin3L(:,1), Bilin3L(:,2),'LineWidth', 3, 'color', 'yellow');
ylabel('[m^3]', 'FontSize', 28);
xlabel('Time [days]', 'FontSize', 28);
title('Oil cumulative production [m^3]', 'FontSize', 28);
legend({'Fine scale','Const_{2L}', 'Bilin_{2L}','Const_{3L}','Bilin_{3L}'}, 'FontSize', 28, 'FontWeight', 'bold', 'Location', 'best');
set(gca,'fontsize',28)

%% Compute error
ErrorConst2L = abs(Const2L(end, 2) - FineGrid(end,2))/FineGrid(end,2)*100
ErrorBilin2L = abs(Bilin2L(end, 2) - FineGrid(end,2))/FineGrid(end,2)*100
ErrorConst3L = abs(Const3L(end, 2) - FineGrid(end,2))/FineGrid(end,2)*100
ErrorBilin3L = abs(Bilin3L(end, 2) - FineGrid(end,2))/FineGrid(end,2)*100