
function PlotProductionCurves(CaseDirectory, color, well)
Ph1_file = strcat(CaseDirectory, '/Output/Prod_Phase1.txt');
Ph2_file = strcat(CaseDirectory, '/Output/Prod_Phase2.txt');
Comp1_file = strcat(CaseDirectory, '/Output/Prod_Comp1.txt');
Comp2_file = strcat(CaseDirectory, '/Output/Prod_Comp2.txt');

Ph1_prod = load(Ph1_file);
Ph2_prod = load(Ph2_file);
% Comp1_prod = load(Comp1_file);
% Comp2_prod = load(Comp2_file);
well = well + 1;

figure(1)
title('Prod ph1');
plot(Ph1_prod(:,1), Ph1_prod(:,well), color);
legend('FS','ADM','UPS')
hold on

figure(2)
title('Brine Production');
plot(Ph2_prod(:,1), Ph2_prod(:,well), color);
xlabel('Days')
ylabel('PV production')
%legend('FS','ADM','UPS')
[hleg1, hobj1] = legend('FS','ADM','UPS');
textobj = findobj(hobj1, 'type', 'text');
set(textobj, 'Interpreter', 'latex', 'fontsize', 15);
hold on



% figure(3)
% title('Prod comp1');
% plot(Comp1_prod(:,1), Comp1_prod(:,end), color);
% hold on
% figure(4)
% title('Prod comp2');
% plot(Comp2_prod(:,1), Comp2_prod(:,end), color);
% hold on
end