%Plot for 1D problems
    x=linspace(Grid.Lx/(2*Grid.Nx), (2*Grid.Nx*Grid.Lx-Grid.Lx)/(2*Grid.Nx), Grid.Nx);
    figure(4)
    subplot(2,1,1);
    plot(x, P, 'green', 'LineWidth',3);
    hold on
    plot(x, P_cap-Pc_cap, 'blue', 'LineWidth',3);
    hold on
    plot(x, P_cap, 'red', 'LineWidth',3);
    legend('no cap','water','oil', 'Location', 'east');
    %title('Pressure [Pa]');
    xlabel('x [m]');
    ylabel('Pressure [Pa]');
    axis([0 Grid.Lx min(P_cap-Pc_cap)-100 max(P_cap)+100])
    
    set(gca,'fontsize',24);
    hold on
    subplot(2,1,2);
    plot(x, 1-S, 'red','LineWidth',3);
    hold on
    plot(x, 1-S_cap, 'blue','LineWidth',3);
    axis([0 Grid.Lx 0 1]);
    xlabel('x [m]');
    ylabel('Oil sat');
    legend('no cap','cap', 'Location', 'east');
    set(gca,'fontsize',24);
    drawnow;