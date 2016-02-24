function PcPlot(Fluid)
S=0:0.01:1;
[Pc, dPc] = ComputePc(S, Fluid);
figure(125);
subplot(1,2,1);
plot(S, Pc);
xlabel('water saturation');
ylabel('Pcow [Pa]');
subplot(1,2,2);
plot(S, dPc);
xlabel('water saturation');
ylabel('dPc/dS');
end