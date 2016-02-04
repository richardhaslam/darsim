% Scalability analysis of ADM for homogeneous problem
Problem = 'Homogeneous';
%% Read data from files
%Stats
file = strcat('../ResultsAndImages/Scalability/',Problem,'/HomoSmall/DLGRStat.txt');
HomoSmall.Stats = load(file);
file = strcat('../ResultsAndImages/Scalability/',Problem,'/HomoMedium/DLGRStat.txt');
HomoMedium.Stats = load(file);
file = strcat('../ResultsAndImages/Scalability/',Problem,'/HomoLarge/DLGRStat.txt');
HomoLarge.Stats = load(file);

%% Active nodes: print Data file for pgfplots
Small_fine = 108 * 27;
Medium_fine = 216 * 54;
Large_fine = 432 * 108;
HomoSmall_Nodes = mean(sum(HomoSmall.Stats(:,4:6),2));
HomoMedium_Nodes = mean(sum(HomoMedium.Stats(:,4:6),2));
HomoLarge_Nodes = mean(sum(HomoLarge.Stats(:,4:6),2));
HomoSmall_AvNodes = mean(sum(HomoSmall.Stats(:,4:6),2)./Small_fine*100);
HomoMedium_AvNodes = mean(sum(HomoMedium.Stats(:,4:6),2)./Medium_fine*100);
HomoLarge_AvNodes = mean(sum(HomoLarge.Stats(:,4:6),2)./Large_fine*100);
%file = fopen(strcat('../ResultsAndImages/Scalability/',Problem,'/ADM/HomoSmall/HomoSmall.dat'), 'w');
%fprintf(file, '%10.5f %10.5f\n', [HomoSmall.Stats(:,1)/HomoSmall(End,1), HomoSmall.Stats(:,4)]');
%fclose(file);
%file = fopen(strcat('../ResultsAndImages/Scalability/',Problem,'/ADM/HomoMedium/HomoMedium.dat'), 'w');
%fprintf(file, '%10.5f %10.5f\n', [HomoMedium.Stats(:,1)/HomoMedium(End,1), HomoMedium.Stats(:,4)]');
%fclose(file);
%file = fopen(strcat('../ResultsAndImages/Scalability/',Problem,'/ADM/HomoLarge/HomoLarge.dat'), 'w');
%fprintf(file, '%10.5f %10.5f\n', [HomoSmall.Stats(:,1)/HomoSmall(End,1), HomoSmall.Stats(:,4)]');
%fclose(file);

%% Number of nonlinear iterations
HomoSmall_nl = mean(HomoSmall.Stats(:,3));
HomoMedium_nl = mean(HomoMedium.Stats(:,3));
HomoLarge_nl = mean(HomoLarge.Stats(:,3));
% file = fopen(strcat('../ResultsAndImages/Scalability/',Problem,'/ADM/', DS, '/',Problem,'_', DS, '.dat'), 'w');
% fprintf(file, '%20s %10s %10s %10s \n','', 'HomoSmall', 'HomoMedium', 'HomoLarge');
% fprintf(file, '%20s %10.2f %10.2f %10.2f \n', '% active nodes', [mean(sum(HomoSmall.Stats(:,4:6), 2)./FineCells(1)*100),...
%     mean(sum(HomoMedium.Stats(:,4:6), 2)./FineCells(2)*100), mean(sum(HomoLarge.Stats(:,4:6),2)./FineCells(3)*100)]);
% % fprintf(file, '%20s %10.0f %10.0f %10.0f \n', 'NnLinear', [HomoSmall_nl, HomoMedium_nl, HomoLarge_nl]);
%fprintf(file, '%20s %10.5f %10.5f %10.5f \n', 'Pressure error', [HomoSmall_EP, HomoMedium_EP, HomoLarge_EP]);
%fprintf(file, '%20s %10.5f %10.5f %10.5f \n', 'Saturation error', [HomoSmall_ES, HomoMedium_ES, HomoLarge_ES]);
%fclose(file);

%% CPU Time
file = strcat('../ResultsAndImages/Scalability/',Problem,'/HomoSmall/DLGRTimings.txt');
HomoSmall.Timings = load(file);
file = strcat('../ResultsAndImages/Scalability/',Problem,'/HomoMedium/DLGRTimings.txt');
HomoMedium.Timings = load(file);
file = strcat('../ResultsAndImages/Scalability/',Problem,'/HomoLarge/DLGRTimings.txt');
HomoLarge.Timings = load(file);
file = strcat('../ResultsAndImages/Scalability/',Problem,'/HomoSmall/FIMTimings.txt');
HomoSmallFine.Timings = load(file);
file = strcat('../ResultsAndImages/Scalability/',Problem,'/HomoMedium/FIMTimings.txt');
HomoMediumFine.Timings = load(file);
file = strcat('../ResultsAndImages/Scalability/',Problem,'/HomoLarge/FIMTimings.txt');
HomoLargeFine.Timings = load(file);

%Linear Solver
HomoSmall_LS = mean(HomoSmall.Timings(:,5))/mean(HomoSmall.Timings(:,5));
HomoMedium_LS = mean(HomoMedium.Timings(:,5))/mean(HomoSmall.Timings(:,5));
HomoLarge_LS = mean(HomoLarge.Timings(:,5))/mean(HomoSmall.Timings(:,5));
HomoSmall_fineLS = mean(HomoSmallFine.Timings(:,4))/mean(HomoSmallFine.Timings(:,4));
HomoMedium_fineLS = mean(HomoMediumFine.Timings(:,4))/mean(HomoSmallFine.Timings(:,4));
HomoLarge_fineLS = mean(HomoLargeFine.Timings(:,4))/mean(HomoSmallFine.Timings(:,4));

HomoSmall_tot = mean(HomoSmall.Timings(:,2))/mean(HomoSmall.Timings(:,2))
HomoMedium_tot = mean(HomoMedium.Timings(:,2))/mean(HomoSmall.Timings(:,2))
HomoLarge_tot = mean(HomoLarge.Timings(:,2))/mean(HomoSmall.Timings(:,2))
HomoSmall_finetot = mean(HomoSmallFine.Timings(:,2))/mean(HomoSmall.Timings(:,2))
HomoMedium_finetot = mean(HomoMediumFine.Timings(:,2))/mean(HomoSmall.Timings(:,2))
HomoLarge_finetot = mean(HomoLargeFine.Timings(:,2))/mean(HomoSmall.Timings(:,2))

%
HomoSmall_RPconstruct = mean(HomoSmall.Timings(:,3))/mean(HomoSmall.Timings(:,3));
HomoMedium_RPconstruct = mean(HomoMedium.Timings(:,3))/mean(HomoSmall.Timings(:,3));
HomoLarge_RPConstruct = mean(HomoLarge.Timings(:,3))/mean(HomoSmall.Timings(:,3));

