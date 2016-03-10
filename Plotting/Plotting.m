%Plotting - Pressure and Saturation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pressure_3D=Options.Pressure_3D;
problem_1D=Options.problem_1D;
if (problem_1D==1)
    %Plot for 1D problems
    x=linspace(Grid.Lx/(2*Grid.Nx), (2*Grid.Nx*Grid.Lx-Grid.Lx)/(2*Grid.Nx), Grid.Nx);
    figure(4)
    subplot(2,1,1);
    plot(x, P, 'red', 'LineWidth',5);
    hold on
    plot(x, P-Pc, 'blue', 'LineWidth',5);
    %title('Pressure [Pa]');
    xlabel('x [m]');
    ylabel('Pressure [Pa]');
    axis([0 Grid.Lx min(P-Pc)-100 max(P)+100])
    legend('oil','water', 'Location', 'east');
     set(gca,'fontsize',24);
    hold on
    subplot(2,1,2);
    plot(x, S, 'red', 'LineWidth', 3);
    axis([0 Grid.Lx 0 1]);
    title('Saturation of water');
    xlabel('x [m]');
    ylabel('Saturation');
    hold on;
    set(gca,'fontsize',24);
    drawnow;
else
    %Plot for 2D problems
    x=linspace(Grid.Lx/(2*Grid.Nx), (2*Grid.Nx*Grid.Lx-Grid.Lx)/(2*Grid.Nx), Grid.Nx);
    y=linspace(Grid.Ly/(2*Grid.Ny), (2*Grid.Ny*Grid.Ly-Grid.Ly)/(2*Grid.Ny), Grid.Ny);
    [X, Y] = meshgrid(x,y);
    %Contour Plot
    if (Options.ContourPlot==1)
        %figure(4)
        %suptitle([num2str(round(t/(3600*24))) ' days']);
        
        %Pressure plot
        %subplot(2,1,1);
        figure(300)
        if (Pressure_3D==1)
            surf(X,Y,P,P);
        else
            contourf(x,y,P');
            if (Grid.Nx==Grid.Ny)
                axis square;
            end
        end
        colorbar;
        xlabel('x [m]');
        ylabel('y [m]');
        c=colorbar;
        pressuremax =max(max(P));
        caxis ([0.99*pressuremax,pressuremax]);
        axis ([0 100 0 100 0.99*pressuremax pressuremax]);
        ylabel(c,'p [Pa]');
        view(75, 35);
        
        %Saturation plot
        %subplot(2,1,2);
        figure(400)
        surf(X,Y,S,S);
        if (Grid.Nx==Grid.Ny)
            axis square;
        end
        caxis ([0 1]);
        colorbar;
        title('Saturation of water');
        xlabel('x [m]');
        ylabel('y [m]');
        view(75, 35);
        drawnow
    else
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
        pressuremax =max(max(P));
        pressuremin = min(min(P));
        caxis ([Prod(1).p,Inj(1).p]);
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
        drawnow
        
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
        caxis ([pressuremin, pressuremax]);
        axis('image');
        set(gca,'fontsize',24);
        
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