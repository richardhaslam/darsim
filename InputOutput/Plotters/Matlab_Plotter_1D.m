% Matlab Plotter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 11 July 2016
%Last modified: 11 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Matlab_Plotter_1D < Plotter
    properties
        color1
        color2
    end
    methods
        function obj = Matlab_Plotter_1D()
            obj.color1 = 'red';
            obj.color2 = 'green';
            obj.VTKindex = 1;
        end
        function PlotSolution(obj, ProductionSystem, DiscretizationModel)
            Grid = DiscretizationModel.ReservoirGrid;
            Status = ProductionSystem.Reservoir.State;
            %Plot for 1D problems
            x = linspace(Grid.Nx * Grid.dx/(2*Grid.Nx), (2*Grid.Nx^2*Grid.dx-Grid.Nx * Grid.dx)/(2*Grid.Nx), Grid.Nx);
            
            N_var = double(Status.Properties.Count);
            Names = Status.Properties.keys;
            for i=1:N_var
                if Status.Properties(Names{i}).Plot
                    figure(i+1)
                    plot(x, Status.Properties(Names{i}).Value, obj.color1, 'LineWidth',1);
                    xlabel('x [m]');
                    ylabel(Names{i});
                    %axis([0 Grid.dx*Grid.Nx Status.Properties(Names{i}).Valmin Status.Properties(Names{i}).Valmax])
                    set(gca,'fontsize',24);
                    hold on
                    drawnow;
                end
            end
        end
        function PlotPermeability(obj, Grid, Perm)
            figure(1)
            x = linspace(0, Grid.Nx*Grid.dx-Grid.dx, Grid.Nx);;
            plot(x, log10(Perm));
            %title('Permeability');
            axis([0 Grid.dx*Grid.Nx min(min(log10(Perm)))-1 max(max(log10(Perm)))+1]);
            %axis([0 Grid.Lx 0 1]);
            xlabel('x [m]');
            ylabel('Log(K) [m^2]');
            set(gca,'fontsize',24);
        end
        function PlotResidual(obj)
        end
        function PlotWells(obj, Inj, Prod, Grid)
        end
        function PlotADMGrid(obj, Grid, CoarseGrid)
            
        end
        function PlotDynamicBasisFunctions(obj, ReservoirGrid, ADMProlp)
        end
    end
end