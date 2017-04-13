% ECMOR PAPER: ERRORS vs DS
Directory = '../../ResultsAndImages/ECMORPaper/';
%% Read data from files
%Case 1
file = strcat(Directory,'Case1/MS/DS010/Output/ADMStat.txt');
Case1.Stat = load(file);
file = strcat(Directory,'Case1/MS/DS010/Output/ADMNwProd.txt');
Case1.Prod.Nw = load(file);
file = strcat(Directory,'Case1/MS/DS010/Output/ADMWProd.txt');
Case1.Prod.W = load(file);

%Case 2
file = strcat(Directory,'Case2/MS/DS015/Output/ADMStat.txt');
Case2.Stat = load(file);
file = strcat(Directory,'Case2/MS/DS015/Output/ADMNwProd.txt');
Case2.Prod.Nw = load(file);
file = strcat(Directory,'Case2/MS/DS015/Output/ADMWProd.txt');
Case2.Prod.W = load(file);

%Case 3
file = strcat(Directory,'Case3/MS/DS015/Output/ADMStat.txt');
Case3.Stat = load(file);
file = strcat(Directory,'Case3/MS/DS015/Output/ADMNwProd.txt');
Case3.Prod.Nw = load(file);
file = strcat(Directory,'Case3/MS/DS015/Output/ADMWProd.txt');
Case3.Prod.W = load(file);

Case1.PercNodes = sum(Case1.Stat(:,4:6),2)/(99*99)*100;
Case2.PercNodes = sum(Case2.Stat(:,4:6),2)/(99*99)*100;
Case3.PercNodes = sum(Case3.Stat(:,4:6),2)/(99*99)*100;

Case1.AverageNodes = mean(Case1.PercNodes);
Case2.AverageNodes = mean(Case2.PercNodes);
Case3.AverageNodes = mean(Case3.PercNodes);

Case1.TotalProd = sum(Case1.Prod.Nw(:,2:5), 2) + sum(Case1.Prod.W(:,2:5), 2);
Case2.TotalProd = sum(Case2.Prod.Nw(:,2:5), 2) + sum(Case2.Prod.W(:,2:5), 2);
Case3.TotalProd = sum(Case3.Prod.Nw(:,2:5), 2) + sum(Case3.Prod.W(:,2:5), 2);

%% Print tables for Latex
fileID = fopen(strcat(Directory, 'Case1NodesVsPV.txt'), 'w');
Print = [Case1.TotalProd(2:end),Case1.PercNodes];
fprintf(fileID, '%5.3f %10.6f\n', Print');
fclose(fileID);
fileID = fopen(strcat(Directory, 'Case2NodesVsPV.txt'), 'w');
Print = [Case2.TotalProd(2:end),Case2.PercNodes];
fprintf(fileID, '%5.3f %10.6f\n', Print');
fclose(fileID);
fileID = fopen(strcat(Directory, 'Case3NodesVsPV.txt'), 'w');
Print = [Case3.TotalProd(2:end),Case3.PercNodes];
fprintf(fileID, '%5.3f %10.6f\n', Print');
fclose(fileID);

