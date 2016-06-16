%Plotting - Pressure and Saturation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Created: 2015
%Last modified: 7 April 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Plotting(Grid, Status, Fluid, color1, color2, Inj, Prod)

P = reshape(Status.p,Grid.Nx,Grid.Ny);
S = reshape(Status.s,Grid.Nx,Grid.Ny);
Pc = reshape(Status.pc,Grid.Nx,Grid.Ny);
x1 = reshape(Status.x1(:,1),Grid.Nx,Grid.Ny);
x2 = reshape(1-Status.x1(:,2),Grid.Nx,Grid.Ny);
z = reshape(Status.z(:,1),Grid.Nx,Grid.Ny);

problem_1D = 0;
if (Grid.Nx == 1 || Grid.Ny == 1)
            problem_1D = 1;
end
if problem_1D
    %Plot for 1D problems
    x = linspace(Grid.Lx/(2*Grid.Nx), (2*Grid.Nx*Grid.Lx-Grid.Lx)/(2*Grid.Nx), Grid.Nx);
    figure(4)
    subplot(3,1,2);
    plot(x, P, color1, 'LineWidth',1);
    hold on
    plot(x, P-Pc, color2, 'LineWidth',1);
    %title('Pressure [Pa]');
    xlabel('x [m]');
    ylabel('Pressure [Pa]');
    axis([0 Grid.Lx min(P-Pc)-Inj.p*0.1 max(P)+Inj.p*0.1])
    %legend('oil','water', 'Location', 'east');
     set(gca,'fontsize',24);
    hold on
    subplot(3,1,3);
    plot(x, S, color2, 'LineWidth', 1);
    axis([0 Grid.Lx 0 1.1]);
    %title('Saturation of water');
    xlabel('x [m]');
    ylabel('Saturation');
    hold on;
    set(gca,'fontsize',24);
    drawnow;
    
    figure(5)
    subplot(3,1,1);
    plot(x, x1, color2, 'LineWidth', 1);
    axis([0 Grid.Lx .5 1]);
    %title('Component 1 in Phase 1 [-]');
    xlabel('x [m]');
    ylabel('x1w [-]');
    hold on;
    set(gca,'fontsize',24);
    subplot(3,1,2);
    plot(x, x2, color2, 'LineWidth', 1);
    axis([0 Grid.Lx .5 1]);
    %title('Component 2 in Phase 2 [-]');
    xlabel('x [m]');
    ylabel('x2nw [-]');
    hold on;
    set(gca,'fontsize',24);
    subplot(3,1,3);
    plot(x, z, color2, 'LineWidth', 1);
    axis([0 Grid.Lx 0 1]);
    %title('Component 1 Total Mass Fraction [-]');
    xlabel('x [m]');
    ylabel('z1 [-]');
    hold on;
    set(gca,'fontsize',24);
    drawnow;
else
    %Plot for 2D problems
    x = linspace(Grid.Lx/(2*Grid.Nx), (2*Grid.Nx*Grid.Lx-Grid.Lx)/(2*Grid.Nx), Grid.Nx);
    y = linspace(Grid.Ly/(2*Grid.Ny), (2*Grid.Ny*Grid.Ly-Grid.Ly)/(2*Grid.Ny), Grid.Ny);
    [X, Y] = meshgrid(x,y);
    % Grid plot
    %Pressure
    figure(300)
    %suptitle([num2str(round(t/(3600*24))) ' days']);
    
    %subplot(2,1,1);
    h = pcolor(X,Y,P');
    set(h, 'EdgeColor', 'none');
    %view([45 45]);
    if Grid.Nx==Grid.Ny
        axis square;
    end
    title('Pressure [Pa]');
    xlabel('x [m]');
    ylabel('y [m]');
    colormap(jet);
    colorbar;
    [Prod(1).p(1),Inj(1).p(1)]
    caxis ([Prod(1).p(1),Inj(1).p(1)]);
    axis('image');
    set(gca,'fontsize',24);
    
    %Saturation
    figure(400)
    %subplot(2,1,2);
    h = pcolor(X,Y,S');
    set(h, 'EdgeColor', 'none');
    %view([45 45]);
    if Grid.Nx==Grid.Ny
        axis square;
    end
    if (strcmp(Fluid.RelPerm, 'Foam')==1)
        title('Gas Saturation');
    else
        title('Water Saturation');
    end
    xlabel('x [m]');
    ylabel('y [m]');
    colormap(jet);
    colorbar;
    caxis ([0,1]);
    axis('image');
    set(gca,'fontsize',24);
    
    %Capillary pressure
    figure(500)
    h = pcolor(X,Y,Pc');
    set(h, 'EdgeColor', 'none');
    %view([45 45]);
    if Grid.Nx==Grid.Ny
        axis square;
    end
    title('Pc [Pa]');
    xlabel('x [m]');
    ylabel('y [m]');
    colormap(jet);
    colorbar;
    pressuremax =max(max(Pc));
    pressuremin = min(min(Pc));
    %caxis ([pressuremin, pressuremax]);
    axis('image');
    set(gca,'fontsize',24);
    
    drawnow
    %Apparent viscosity
    if (strcmp(Fluid.RelPerm, 'Foam')==1)
        figure(500)
        %subplot(2,1,2)
        [~, ~, ~, ~, app]=Mobilities(S, Fluid);
        h = pcolor(X,Y,app');
        set(h, 'EdgeColor', 'none');
        %view([45 45]);
        if Grid.Nx==Grid.Ny
            axis square;
        end
        title('Apparent Viscosity');
        xlabel('x [m]');
        ylabel('y [m]');
        colormap(jet);
        colorbar;
        caxis ([0,0.09]);
        %caxis ([0,1]);
        axis('image');
        set(gca,'fontsize',24);
        drawnow
    end
end
end
