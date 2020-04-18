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
        function DefinePerforatedCells(obj, Wells)
            % Injectors
            Well_Type = 'injector';
            for w = 1:Wells.NofInj
                obj.DefinePerforatedCellsForEachWell(Wells.Inj(w), Well_Type, w);
            end
            % Producers
            Well_Type = 'producer';
            for w = 1:Wells.NofProd
                obj.DefinePerforatedCellsForEachWell(Wells.Prod(w), Well_Type, w);
            end
        end
        function DefinePerforatedCellsForEachWell(obj, Well, Well_Type, w)
            switch Well.Coordinate.Type
                case('IJK_LIST')
                    obj.ObtainPerforatedCellsBasedOnIJKList(Well, Well_Type, w);
                case('XYZ_LIST')
                    obj.ObtainPerforatedCellsBasedOnXYZList(Well, Well_Type, w);
                case('CELL_INDEX_LIST')
                    obj.ObtainPerforatedCellsBasedOnCellIndexList(Well, Well_Type, w);
                otherwise
                    error('DARSim Error: In of the %s #%d, the coordination keyword "IJK_LIST" or "XYZ_LIST" or "CELL_INDEX_LIST" is missing. Check the input file!\n', Well_Type, w);
            end
            if isempty(Well.Cells)
                error('DARSim Error: The %s #%d does not perforate any grid cell. Please check coordinate of this well in the input file.', Well_Type, w);
            end
            Well.PI = Well.PI(Well.Cells);
        end
        function AverageMassOnCoarseBlocks(obj, Status, FluidModel, Formulation)
            % virtual call
        end
    end
end