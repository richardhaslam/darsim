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
        Pmax
        Pmin
    end
    methods
        function obj = Matlab_Plotter_1D(pmin, pmax)
            obj.color1 = 'green';
            obj.color2 = 'green';
           
            obj.Pmax = pmax*1.1;
            obj.Pmin = pmin*0.9;
        end
        function PlotSolution(obj, Status, Grid)
            % Reshape all objects
            P = reshape(Status.p,Grid.Nx,Grid.Ny);
            S = reshape(Status.S,Grid.Nx,Grid.Ny);
            Pc = reshape(Status.pc,Grid.Nx,Grid.Ny);
            x1 = reshape(Status.x(:,1),Grid.Nx,Grid.Ny);
            x2 = reshape(1-Status.x(:,2),Grid.Nx,Grid.Ny);
            z = reshape(Status.z(:,1),Grid.Nx,Grid.Ny);
            %Plot for 1D problems
            x = linspace(Grid.Nx * Grid.dx/(2*Grid.Nx), (2*Grid.Nx^2*Grid.dx-Grid.Nx * Grid.dx)/(2*Grid.Nx), Grid.Nx);
            figure(2)
            %subplot(3,1,2);
            plot(x, P, obj.color1, 'LineWidth',1);
            hold on
            plot(x, P-Pc, obj.color2, 'LineWidth',1);
            %title('Pressure [Pa]');
            xlabel('x [m]');
            ylabel('Pressure [Pa]');
            axis([0 Grid.dx*Grid.Nx obj.Pmin obj.Pmax])
            %legend('oil','water', 'Location', 'east');
            set(gca,'fontsize',24);
            hold on
            figure(3)
            %subplot(3,1,3);
            plot(x, S, obj.color2, 'LineWidth', 1);
            axis([0 Grid.Nx*Grid.dx 0 1.1]);
            %title('Saturation of water');
            xlabel('x [m]');
            ylabel('Saturation');
            hold on;
            set(gca,'fontsize',24);
            drawnow;
            
            figure(4)
            %subplot(3,1,1);
            plot(x, x1, '.', 'Color', obj.color2);
            axis([0 Grid.Nx*Grid.dx 0 1]);
            %title('Component 1 in Phase 1 [-]');
            xlabel('x [m]');
            ylabel('x1v [-]');
            hold on;
            set(gca,'fontsize',24);
            figure(5)
            %subplot(3,1,2);
            plot(x, x2, '.', 'Color', obj.color2);
            axis([0 Grid.Nx*Grid.dx 0 1]);
            %title('Component 2 in Phase 2 [-]');
            xlabel('x [m]');
            ylabel('x2l [-]');
            hold on;
            set(gca,'fontsize',24);
            %subplot(3,1,3);
            figure(6)
            plot(x, z, obj.color2, 'LineWidth', 1);
            axis([0 Grid.Nx*Grid.dx 0 1]);
            %title('Component 1 Total Mass Fraction [-]');
            xlabel('x [m]');
            ylabel('z1 [-]');
            hold on;
            set(gca,'fontsize',24);
            drawnow;
            obj.color1 = 'red';
            obj.color2 = 'red';
        end
        function PlotPermeability(obj, Grid, Perm)
            figure(1)
            x = linspace(0, Grid.Nx*Grid.dx-Grid.dx, Grid.Nx);
            subplot(3,1,1);
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
    end
end