%Plot capillary pressure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Created: 22 April 2016
%Last modified: 22 April 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotPc(Pc, dPc, S)
subplot(1,2,1);
plot(S, Pc, 'red', 'LineWidth',2);
xlabel('water saturation');
ylabel('Pcow [Pa]');
axis('square');
set(gca,'fontsize',24);
hold on
subplot(1,2,2);
plot(S, dPc, 'red', 'LineWidth',2);
xlabel('water saturation');
ylabel('dPc/dS');
axis('square');
set(gca,'fontsize',24);
hold on
end