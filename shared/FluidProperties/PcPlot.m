function PcPlot(Fluid)
S=0:0.01:1;
[Pc, dPc] = ComputePc(S, Fluid);
figure(125);
%subplot(1,2,1);
plot(S, Pc, 'LineWidth',5);
xlabel('water saturation');
ylabel('Pcow [Pa]');
axis('square');
set(gca,'fontsize',24);
%subplot(1,2,2);
figure(126)
plot(S, dPc, 'LineWidth',5);
xlabel('water saturation');
ylabel('dPc/dS');
axis('square');
set(gca,'fontsize',24);
end