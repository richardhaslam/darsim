% Reservoir Simulator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
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
            obj.Builder = simulation_builder();
            obj.Simulation = Reservoir_Simulation();
        end
        function BuildObjects(obj)
            obj.Simulation = obj.Builder.BuildSimulation(obj.Reader.FractureMatrix);
            obj.Writer = obj.Builder.BuildWriter(obj.Reader.Directory, obj.Simulation); 
        end
        function PrintInfo(obj)
            disp(['******************', num2str(obj.Builder.SimulationInput.ProblemName),'******************']);
            disp(newline);
            disp('FORMATION CHARACTERISTICS:');
           
            disp('Reservoir geometry:');
            disp(['Lx: ', num2str(obj.Simulation.ProductionSystem.Reservoir.Length), ' m']);
            disp(['Ly: ', num2str(obj.Simulation.ProductionSystem.Reservoir.Width), ' m']);
            disp(['Thickness:  ', num2str(obj.Simulation.ProductionSystem.Reservoir.Thickness), ' m']);
            Nx = obj.Simulation.DiscretizationModel.ReservoirGrid.Nx;
            Ny = obj.Simulation.DiscretizationModel.ReservoirGrid.Ny;
            Nz = obj.Simulation.DiscretizationModel.ReservoirGrid.Nz;
            disp(['Total Grid Cells: ', num2str(Nx), ' x ',  num2str(Ny), ' x ', num2str(Nz), ' = ', num2str(Nx*Ny*Nz), ', Number of Active Cells: ',num2str(obj.Simulation.DiscretizationModel.ReservoirGrid.N)]);
            disp('---------------------------------------------------------');
            if obj.Builder.SimulationInput.FracturesProperties.isFractured
                for f = 1:obj.Simulation.DiscretizationModel.FracturesGrid.Nfrac
                    Nx = obj.Simulation.DiscretizationModel.FracturesGrid.Grids(f).Nx;
                    Ny = obj.Simulation.DiscretizationModel.FracturesGrid.Grids(f).Ny;
                    fprintf('Fracture %2d: Grid= %3.0f x %3.0f = %3.0f\n', f, Nx, Ny, Nx*Ny);
                end
                fprintf('Total fracture grids: %3.0f\n', sum(obj.Simulation.DiscretizationModel.FracturesGrid.N));
                disp('---------------------------------------------------------');
            end
            disp('Fluid Model:');
            disp(['Type: ', obj.Simulation.FluidModel.name]);
            disp(['N of phases: ', num2str(obj.Simulation.FluidModel.NofPhases)]);
            disp(['N of components:  ', num2str(obj.Simulation.FluidModel.NofComponents)]);
            disp('---------------------------------------------------------');
            
            Time = obj.Simulation.TimeDriver.Sec2DHMS(obj.Simulation.TimeDriver.TotalTime);
            disp(['Total Simulation Time: ' num2str(obj.Simulation.TimeDriver.TotalTime) ' seconds (', ...
                num2str(Time.Days), ' days : ', num2str(Time.Hours), ' hrs : ', num2str(Time.Minutes), ' mins : ', num2str(Time.Seconds), ' sec)']);
            disp(['Maximum Allowed Time-steps: ', num2str(obj.Builder.SimulatorSettings.MaxNumTimeSteps) ]);
            disp(['Minimum Size of Time-step : ', num2str(obj.Builder.SimulatorSettings.MinMaxdt(1)) ]);
            disp(['Maximum Size of Time-step : ', num2str(obj.Builder.SimulatorSettings.MinMaxdt(2)) ]);
            disp(['Number of Reports         : ', num2str(obj.Builder.SimulatorSettings.Reports)     ]);
            disp('---------------------------------------------------------');
        end
        function Run(obj)
            % Initialize Simulation
            obj.Simulation.Initialize();
            % Plot initial state of the reservoir
            obj.Writer.PlotSolution(obj.Simulation.ProductionSystem, obj.Simulation.DiscretizationModel);
            % Write initial state on a file
            obj.Writer.WriteSolutionOnFile(obj.Simulation.ProductionSystem)
            obj.Writer.Index = obj.Writer.Index + 1;
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