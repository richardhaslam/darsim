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
            obj.color1 = 'green';
            obj.color2 = 'green';
        end
        function PlotSolution(obj, Status, Grid)
            % Reshape all objects
            P = reshape(Status.p,Grid.Nx,Grid.Ny);
            S = reshape(Status.s,Grid.Nx,Grid.Ny);
            Pc = reshape(Status.pc,Grid.Nx,Grid.Ny);
            x1 = reshape(Status.x1(:,1),Grid.Nx,Grid.Ny);
            x2 = reshape(1-Status.x1(:,2),Grid.Nx,Grid.Ny);
            z = reshape(Status.z(:,1),Grid.Nx,Grid.Ny);
            %Plot for 1D problems
            x = linspace(Grid.Lx/(2*Grid.Nx), (2*Grid.Nx*Grid.Lx-Grid.Lx)/(2*Grid.Nx), Grid.Nx);
            figure(1)
            subplot(3,1,2);
            plot(x, P, obj.color1, 'LineWidth',1);
            hold on
            plot(x, P-Pc, obj.color2, 'LineWidth',1);
            %title('Pressure [Pa]');
            xlabel('x [m]');
            ylabel('Pressure [Pa]');
            axis([0 Grid.Lx min(P-Pc)-max(P)*0.1 max(P)+max(P)*0.1])
            %legend('oil','water', 'Location', 'east');
            set(gca,'fontsize',24);
            hold on
            subplot(3,1,3);
            plot(x, S, obj.color2, 'LineWidth', 1);
            axis([0 Grid.Lx 0 1.1]);
            %title('Saturation of water');
            xlabel('x [m]');
            ylabel('Saturation');
            hold on;
            set(gca,'fontsize',24);
            drawnow;
            
            figure(2)
            subplot(3,1,1);
            plot(x, x1, obj.color2, 'LineWidth', 1);
            axis([0 Grid.Lx .5 1]);
            %title('Component 1 in Phase 1 [-]');
            xlabel('x [m]');
            ylabel('x1w [-]');
            hold on;
            set(gca,'fontsize',24);
            subplot(3,1,2);
            plot(x, x2, obj.color2, 'LineWidth', 1);
            axis([0 Grid.Lx .5 1]);
            %title('Component 2 in Phase 2 [-]');
            xlabel('x [m]');
            ylabel('x2nw [-]');
            hold on;
            set(gca,'fontsize',24);
            subplot(3,1,3);
            plot(x, z, obj.color2, 'LineWidth', 1);
            axis([0 Grid.Lx 0 1]);
            %title('Component 1 Total Mass Fraction [-]');
            xlabel('x [m]');
            ylabel('z1 [-]');
            hold on;
            set(gca,'fontsize',24);
            drawnow;
            obj.color1 = 'red';
            obj.color2 = 'blue';
        end
        function PlotPermeability(obj, Grid, Perm)
            figure(1)
            x = linspace(0, Grid.Lx-Grid.dx, Grid.Nx);
            subplot(3,1,1);
            plot(x, log10(Perm));
            %title('Permeability');
            axis([0 Grid.Lx min(min(log10(Perm)))-1 max(max(log10(Perm)))+1]);
            %axis([0 Grid.Lx 0 1]);
            xlabel('x [m]');
            ylabel('Log(K) [m^2]');
            set(gca,'fontsize',24);
        end
        function PlotResidual(obj)
        end
        function PlotWells(obj, Inj, Prod)
        end
        function PlotADMGrid(obj, Grid, CoarseGrid)
            
        end
    end
end