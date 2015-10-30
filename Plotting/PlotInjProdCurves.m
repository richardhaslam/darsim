%Injection and Production data

figure (6)
subplot(2,1,1);
plot(CumulativeTime(1:Ndt), Inj.water(1:Ndt));
ylabel('[m^3]');
xlabel('Time [days]')
title('Water cumulative injection [m^3]')

subplot(2,1,2);
%plot(CumulativeTime(1:Ndt), Prod.water(1:Ndt));
hold on
plot(CumulativeTime(1:Ndt), Prod.oil(1:Ndt),'LineWidth', 3, 'color', 'black');
ylabel('[m^3]');
xlabel('Time [days]')
title('Oil cumulative production [m^3]');
legend({'Bilin_{2L}', 'Const_{2L}','BaseGrid'}, 'FontSize', 16, 'FontWeight', 'bold');
