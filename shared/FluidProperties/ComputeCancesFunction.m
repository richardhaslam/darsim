%Compute Cances function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Fluid] = ComputeCancesFunction(Fluid, K, por)
s = 0.1:0.001:1;
%Numerical 
[Mw, Mo] = Mobilities(s, Fluid);
Mt = Mw + Mo;
[~, dPc, dJ] = ComputePc(s, Fluid, K, por);
PcPlot(Pc, dPc, s);
switch (Fluid.Pc)
    case('BrooksCorey')
        integrand = (Mw.*Mo)./Mt.*dPc*(-1);
        C = zeros(length(s)-1, 1);
        for i=2:length(s)
            C(i-1) = - trapz(s(1:i), integrand(1:i));
        end
        Cfitted = fit(s(2:end)',C,'poly5');
        figure(177)
        C = [0;C];
        plot(s(1:end), C, 'red', 'LineWidth', 5);
        xlabel('water saturation');
        ylabel('Cances function');
        hold on
        C2 = Cfitted(s);
        plot(s, C2, 'blue', 'LineWidth', 3);
        legend('Cances','fitted');
        set(gca,'fontsize',24);
        Fluid.Cances = Cfitted;
    case('JLeverett')
         integrand = (Mw.*Mo)./Mt.*dPc*(-1);
        C = zeros(length(s)-1, 1);
        for i=2:length(s)
            C(i-1) = - trapz(s(1:i), integrand(1:i));
        end
        Cfitted = fit(s(2:end)',C,'poly5');
        figure(177)
        C = [0;C];
        plot(s(1:end), C, 'red', 'LineWidth', 5);
        xlabel('water saturation');
        ylabel('Cances function');
        hold on
        C2 = Cfitted(s);
        plot(s, C2, 'blue', 'LineWidth', 3);
        legend('Cances','fitted');
        set(gca,'fontsize',24);
        Fluid.Cances = Cfitted;  
end
end