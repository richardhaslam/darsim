% Matlab Plotter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Matlab_Plotter_2D < Plotter
    properties (Constant)
        FontSize = 24
        
    end
    properties
    end
    methods
        function PlotSolution(obj, ProductionSystem, DiscretizationModel)
            obj.PlotReservoirSolution(ProductionSystem.Reservoir.State, DiscretizationModel.ReservoirGrid);
            for f = 1 : length(ProductionSystem.FracturesNetwork.Fractures)
                obj.PlotFractureSolution(ProductionSystem.FracturesNetwork.Fractures(f), DiscretizationModel.FracturesGrid.Grids(f), f);
            end
        end
        function PlotReservoirSolution(obj, Status, Grid)
            %Plot for 2D problems
            x = linspace(Grid.Nx * Grid.dx/(2*Grid.Nx), (2*Grid.Nx^2*Grid.dx-Grid.Nx * Grid.dx)/(2*Grid.Nx), Grid.Nx);
            y = linspace(Grid.Ny * Grid.dy/(2*Grid.Ny), (2*Grid.Ny^2*Grid.dy-Grid.Ny * Grid.dy)/(2*Grid.Ny), Grid.Ny);
            [X, Y] = meshgrid(x,y);
            
            N_var = double(Status.Properties.Count);
            Names = Status.Properties.keys;
            for i=1:N_var
                figure(i+1)
                Var = reshape(Status.Properties(Names{i}).Value, Grid.Nx, Grid.Ny);
                h = pcolor(X,Y, Var');
                set(h, 'EdgeColor', 'none');
                caxis ([Status.Properties(Names{i}).Valmin, Status.Properties(Names{i}).Valmax]);
                %view([45 45]);
                if Grid.Nx==Grid.Ny
                    axis square;
                end
                title(Names{i});
                xlabel('x [m]');
                ylabel('y [m]');
                colormap(jet);
                colorbar;
                axis('image');
                set(gca,'fontsize',24);
                drawnow
            end
        end
        function PlotFractureSolution(obj, Fracture, Grid, f)
           % TO BE IMPLEMENTED 
        end
        function PlotPermeability(obj, Grid, Perm)
            %Plot permeability
            x = linspace(0, Grid.Nx*Grid.dx-Grid.dx, Grid.Nx);
            y = linspace(0, Grid.Ny*Grid.dy-Grid.dy, Grid.Ny);
            [X, Y] = meshgrid(x,y);
            %Use log scale
            K = log10(Perm(:,1));
            
            figure(1)
            K = reshape(K, Grid.Nx, Grid.Ny);
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
        function PlotWells(obj, Inj, Prod, Grid)
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