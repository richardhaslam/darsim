%  Discretization model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%  Discretization model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Discretization_model < handle
    properties
        N
        ReservoirGrid
        FracturesGrid
        CrossConnections
    end
    methods
        function AddReservoirGrid(obj, reservoirgrid)
            obj.ReservoirGrid = reservoirgrid;
        end
        function AddFracturesGrid(obj, fracturesgrid)
            obj.FracturesGrid = fracturesgrid;
        end
        function AddCrossConnections(obj, crossconnections, Formulation)
            obj.CrossConnections = crossconnections; 
        end
        function Initialize(obj, ProductionSystem, Formulation)
            obj.ReservoirGrid.Initialize(ProductionSystem.Reservoir);
            obj.ReservoirGrid.AddGridCoordinates();
            obj.ReservoirGrid.CorrectTransmissibilitiesForpEDFM();
            % Perforated cells
            obj.DefinePerforatedCells(ProductionSystem.Wells);
            
            % Total number of cells
            obj.N = obj.ReservoirGrid.N;
            
            % . Assign Depth
            obj.ReservoirGrid.ComputeDepth(Formulation.GravityModel.alpha, ProductionSystem.Reservoir.Thickness);
            obj.N = obj.ReservoirGrid.N;
            
            % initialize fractures
            if ProductionSystem.FracturesNetwork.Active
                for f = 1:ProductionSystem.FracturesNetwork.NumOfFrac
                    obj.FracturesGrid.Grids(f).Initialize(ProductionSystem.FracturesNetwork.Fractures(f));
                    obj.FracturesGrid.Grids(f).CorrectTransmissibilitiesForpEDFM();
                end
                % Total number of cells
                obj.N = obj.N + sum(obj.FracturesGrid.N);
                % Adding the harmonic permeabilities to CrossConnections
                obj.AddHarmonicPermeabilities(ProductionSystem.Reservoir, ProductionSystem.FracturesNetwork.Fractures);
                % Adding the harmonic conductivities to CrossConnections
                if ~isempty(ProductionSystem.Reservoir.K_Cond_eff)
                    obj.AddHarmonicConductivities(ProductionSystem.Reservoir, ProductionSystem.FracturesNetwork.Fractures);
                end
            end
        end
        function AverageMassOnCoarseBlocks(obj, Status, FluidModel, Formulation)
            % virtual call
        end
    end
end