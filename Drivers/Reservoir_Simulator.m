% Reservoir Simulator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 12 July 2016
%Last modified: 28 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Reservoir_Simulator < handle
    properties
        Reader
        Builder
        Simulation
        Writer
    end
    methods
        function obj = Reservoir_Simulator(Directory, File)
            obj.Reader = reader_darsim2(Directory, File);
            obj.Builder = builder();
            obj.Simulation = Reservoir_Simulation();
        end
        function BuildObjects(obj)
            obj.Builder.FindKeyWords(obj.Reader.InputMatrix, obj.Reader.SettingsMatrix, obj.Reader.FractureMatrix);
            obj.Simulation = obj.Builder.BuildSimulation(obj.Reader.InputMatrix{1}, obj.Reader.SettingsMatrix{1}, obj.Reader.FractureMatrix);
            obj.Writer = obj.Builder.BuildWriter(obj.Reader.Directory, obj.Simulation); 
        end
        function PrintInfo(obj)
            disp(['******************', num2str(obj.Builder.ProblemName),'******************']);
            disp(newline);
            disp('FORMATION CHARACTERISTICS:');
           
            disp('Reservoir geometry:');
            disp(['Lx: ', num2str(obj.Simulation.ProductionSystem.Reservoir.Length), ' m']);
            disp(['Ly: ', num2str(obj.Simulation.ProductionSystem.Reservoir.Width), ' m']);
            disp(['Thickness:  ', num2str(obj.Simulation.ProductionSystem.Reservoir.Thickness), ' m']);
            disp(['Grid: ', num2str(obj.Simulation.DiscretizationModel.ReservoirGrid.Nx), ' x ',  num2str(obj.Simulation.DiscretizationModel.ReservoirGrid.Ny), ' x ', num2str(obj.Simulation.DiscretizationModel.ReservoirGrid.Nz)]);
            disp('---------------------------------------------------------');
            if obj.Builder.Fractured
                for f = 1:obj.Simulation.DiscretizationModel.FracturesGrid.Nfrac
                    fprintf('Fracture %2d: Grid= %3.0f x %3.0f\n', f, ...
                        obj.Simulation.DiscretizationModel.FracturesGrid.Grids(f).Nx, obj.Simulation.DiscretizationModel.FracturesGrid.Grids(f).Ny)
                end
                fprintf('Total fracture grids: %3.0f\n', sum(obj.Simulation.DiscretizationModel.FracturesGrid.N));
                disp('---------------------------------------------------------');
            end
            disp('Fluid Model:');
            disp(['Type: ', obj.Simulation.FluidModel.name]);
            disp(['N of phases: ', num2str(obj.Simulation.FluidModel.NofPhases)]);
            disp(['N of components:  ', num2str(obj.Simulation.FluidModel.NofComp)]);
            disp('---------------------------------------------------------');
        end
        function Run(obj)
            % Initialize Simulation
            obj.Simulation.Initialize();
            % Plot initial state of the reservoir
            obj.Writer.PlotSolution(obj.Simulation.ProductionSystem, obj.Simulation.DiscretizationModel);
            % Write initial state on a file
            obj.Writer.WriteSolutionOnFile(obj.Simulation.ProductionSystem, 0)
            % Plot Wells
            obj.Writer.Plotter.PlotWells(obj.Simulation.ProductionSystem.Wells.Inj, obj.Simulation.ProductionSystem.Wells.Prod, obj.Simulation.DiscretizationModel.ReservoirGrid);
            % Run simulation
            obj.Simulation.Run(obj.Writer);
        end
        function OutputResults(obj)
            obj.Writer.WriteSummary(obj.Simulation.Summary);
        end
    end
end