% Matlab Plotter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 11 July 2016
%Last modified: 11 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Matlab_Plotter_2D < Plotter
    properties (Constant)
        FontSize = 24
    end
    methods
        function PlotSolution(obj, Status, Grid)
            %Reshape objects for pcolor plots
            P = reshape(Status.p,Grid.Nx,Grid.Ny);
            S = reshape(Status.S,Grid.Nx,Grid.Ny);
            Pc = reshape(Status.pc,Grid.Nx,Grid.Ny);
            %x1 = reshape(Status.x1(:,1),Grid.Nx,Grid.Ny);
            %x2 = reshape(1-Status.x1(:,2),Grid.Nx,Grid.Ny);
            %z = reshape(Status.z(:,1),Grid.Nx,Grid.Ny);
            
            %Plot for 2D problems
            x = linspace(Grid.Nx * Grid.dx/(2*Grid.Nx), (2*Grid.Nx^2*Grid.dx-Grid.Nx * Grid.dx)/(2*Grid.Nx), Grid.Nx);
            y = linspace(Grid.Ny * Grid.dy/(2*Grid.Ny), (2*Grid.Ny^2*Grid.dy-Grid.Ny * Grid.dy)/(2*Grid.Ny), Grid.Ny);
            [X, Y] = meshgrid(x,y);
            % Grid plot
            %Pressure
            figure(1)
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
            axis('image');
            set(gca,'fontsize',24);
            
            %Saturation
            figure(2)
            h = pcolor(X,Y,S');
            set(h, 'EdgeColor', 'none');
            if Grid.Nx==Grid.Ny
                axis square;
            end
          
            xlabel('x [m]');
            ylabel('y [m]');
            colormap(jet);
            colorbar;
            caxis ([0,1]);
            axis('image');
            set(gca,'fontsize', obj.FontSize);
            
            %Capillary pressure
            figure(3)
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
            set(gca,'fontsize', obj.FontSize);
            drawnow
        end
        function PlotPermeability(obj, Grid, Perm)
            %Plot permeability
            x = linspace(0, Grid.Nx*Grid.dx-Grid.dx, Grid.Nx);
            y = linspace(0, Grid.Ny*Grid.dy-Grid.dy, Grid.Ny);
            [X, Y] = meshgrid(x,y);
            %Use log scale
            K(:,:) = Perm(1,:,:);
            K = log10(K);
            
            figure(4)
            K = reshape(Perm(:,1), Grid.Nx, Grid.Ny);
            h = pcolor(X,Y,K');
            set(h, 'EdgeColor', 'none');
            title('Log(K)');
            xlabel('x [m]');
            ylabel('y [m]');
            axis('image');
            colorbar;
            colormap(jet);
            set(gca,'fontsize', obj.FontSize);
            hold on
            drawnow;
        end
        function PlotResidual(obj)
        end
        function PlotWells(obj, Inj, Prod)
        end
        function PlotADMGrid(obj, Grid, CoarseGrid)
            % It s only plotted on top of the Saturation solution
            figure (2)
            hold on
            %Fine Grids
            Nf = Grid.Nx*Grid.Ny;
            for f=1:Nf
                if (Grid.Active(f) == 1)
                    Xmax = Grid.dx*(Grid.I(f));
                    Ymax = Grid.dy*(Grid.J(f));
                    Xmin = Xmax - Grid.dx;
                    Ymin = Ymax - Grid.dy;
                    
                    line([Xmin Xmax], [Ymin Ymin], 'Color', 'w', 'LineWidth',1);
                    line([Xmin Xmax], [Ymax Ymax], 'Color', 'w', 'LineWidth',1);
                    line([Xmin Xmin], [Ymin Ymax], 'Color', 'w', 'LineWidth',1);
                    line([Xmax Xmax], [Ymin Ymax], 'Color', 'w', 'LineWidth',1);
                end
            end
            %Coarse Grid
            for x=1:maxLevel
                Nc = CoarseGrid(x).Nx*CoarseGrid(x).Ny;
                for c=1:Nc
                    if CoarseGrid(x).Active(c) == 1
                        Xmin = Grid.dx*(CoarseGrid(x).I(c) - floor((CoarseGrid(x).CoarseFactor(1) - 1)/2))-Grid.dx;
                        Xmax = Grid.dx*(CoarseGrid(x).I(c) + ceil((CoarseGrid(x).CoarseFactor(1) - 1)/2));
                        Ymin = Grid.dy*(CoarseGrid(x).J(c) - floor((CoarseGrid(x).CoarseFactor(2) - 1)/2))-Grid.dy;
                        Ymax = Grid.dy*(CoarseGrid(x).J(c) + ceil((CoarseGrid(x).CoarseFactor(2) - 1)/2));
                        line([Xmin Xmax], [Ymin Ymin], 'Color', 'w', 'LineWidth',1);
                        line([Xmin Xmax], [Ymax Ymax], 'Color', 'w', 'LineWidth',1);
                        line([Xmin Xmin], [Ymin Ymax], 'Color', 'w', 'LineWidth',1);
                        line([Xmax Xmax], [Ymin Ymax], 'Color', 'w', 'LineWidth',1);
                    end
                end
            end
            drawnow
        end
    end
end