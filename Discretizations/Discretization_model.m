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
        UnifiedGrid
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
                obj.ReservoirGrid.ListOfFracturedReservoirCells = unique(vertcat(obj.CrossConnections.Cells));
                obj.ReservoirGrid.ListOfFracturedReservoirCells(obj.ReservoirGrid.ListOfFracturedReservoirCells>obj.ReservoirGrid.N)=[];
            end
            
            % Adding the non-neighboring connections to the list of neighbors (in case of having fractures)
            if ProductionSystem.FracturesNetwork.Active
                obj.AddNonNeighboringConnectionsToNeighborsListAtFineScale();
            end
        end
        function DefinePerforatedCells(obj, Wells)
            % Injectors
            Well_Type = 'injector';
            for w = 1:Wells.NofInj
                obj.DefinePerforatedCellsForEachWell(Wells.Inj(w), Well_Type, w);
                obj.ReservoirGrid.ListOfPerforatedCells = unique( vertcat(obj.ReservoirGrid.ListOfPerforatedCells,Wells.Inj(w).Cells) );
            end
            % Producers
            Well_Type = 'producer';
            for w = 1:Wells.NofProd
                obj.DefinePerforatedCellsForEachWell(Wells.Prod(w), Well_Type, w);
                obj.ReservoirGrid.ListOfPerforatedCells = unique( vertcat(obj.ReservoirGrid.ListOfPerforatedCells,Wells.Prod(w).Cells) );
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
        function I = Index_Local_to_Global(obj, Im, f, g)
            if (Im>obj.ReservoirGrid.N),  Im = obj.ReservoirGrid.N;  end
            if (f<0),  error('f (fracture index) cannot be negative!');  end
            if (f>obj.FracturesGrid.Nfrac),  error('f exceeds the number of fractures!');  end
            if (g<0),  error('g (fracture cell index) cannot be negative!');  end
            if (f~=0)&&(g>obj.FracturesGrid.Grids(f).N),  error('g exceeds the number of fracture cells!');  end
            
            if (f==0) && (g==0)
                I = Im;
            elseif (f~=0) && (g~=0)
                I = (obj.ReservoirGrid.N) + sum(obj.FracturesGrid.N(1:f-1)) + g;
            else
                error('For global indexing in the reservoir both f & g must be zero!\nFor global indexing in the fractures both f & g must be non-zero!');
            end
        end
        function indexing = Index_Global_to_Local(obj, I)
            if (I<1),  error('Global indexing (I) cannot be negative!');  end
            if (I>obj.N),  error('Global indexing (I) cannot exceed total number of cells!');  end
            if I <= obj.ReservoirGrid.N
                indexing.Im = I;
                indexing.f = 0;
                indexing.g = 0;
            else
                indexing.Im = obj.ReservoirGrid.N;
                temp = I - obj.ReservoirGrid.N;
                temp = find( temp - cumsum(obj.FracturesGrid.N) <= 0);
                indexing.f = temp(1);
                indexing.g = I - obj.ReservoirGrid.N - sum( obj.FracturesGrid.N(1:indexing.f-1) );
                if indexing.g==0,  indexing.g = obj.FracturesGrid.Grids(indexing.f).N;  end
            end
            if obj.Index_Local_to_Global(indexing.Im, indexing.f, indexing.g) ~= I
                error('Im is not correspondent with I. Check the formula again!');
            end
        end
        function AddNonNeighboringConnectionsToNeighborsListAtFineScale(obj)
            % Initializing the NonNeighbours properties as cell arrays
            obj.ReservoirGrid.NonNeighbours = cell(obj.ReservoirGrid.N,1);
            for f = 1 : obj.FracturesGrid.Nfrac
                obj.FracturesGrid.Grids(f).NonNeighbours = cell(obj.FracturesGrid.N(f),1);
            end
            
            % Looping over the connections between each fracture element and its overlaping cells from the other media (reservoir or other fractures)
            for c = 1:length(obj.CrossConnections)
                Frac_Global_Index = c + obj.ReservoirGrid.N;
                Frac_Local_Index = obj.Index_Global_to_Local(Frac_Global_Index);
                f = Frac_Local_Index.f;
                g = Frac_Local_Index.g;
                
                % Assigning matrix-fracture non-neighboring connections
                Reservoir_Overlaps = obj.CrossConnections(c).Cells( obj.CrossConnections(c).Cells <= obj.ReservoirGrid.N);
                obj.FracturesGrid.Grids(f).NonNeighbours{g} = unique([obj.FracturesGrid.Grids(f).NonNeighbours{g}, Reservoir_Overlaps']);
                for n = 1:length(Reservoir_Overlaps)
                    Im = Reservoir_Overlaps(n);
                    obj.ReservoirGrid.NonNeighbours{Im} = unique([obj.ReservoirGrid.NonNeighbours{Im}, Frac_Global_Index]);
                end
                
                % Assigning fracture-fracture non-neighboring connections
                Fracture_Overlaps = obj.CrossConnections(c).Cells( obj.CrossConnections(c).Cells > obj.ReservoirGrid.N);
                obj.FracturesGrid.Grids(f).NonNeighbours{g} = unique([obj.FracturesGrid.Grids(f).NonNeighbours{g}, Fracture_Overlaps']);
                for n = 1:length(Fracture_Overlaps)
                    If = Fracture_Overlaps(n);
                    Frac2_Local_Index = obj.Index_Global_to_Local(If);
                    f2 = Frac2_Local_Index.f;  if f2 == f, error('Sth is wrong'); end
                    g2 = Frac2_Local_Index.g;
                    obj.FracturesGrid.Grids(f2).NonNeighbours{g2} = unique([obj.FracturesGrid.Grids(f2).NonNeighbours{g2}, Frac_Global_Index]);
                end
            end
        end
        function AverageMassOnCoarseBlocks(obj, Status, FluidModel, Formulation)
            % virtual call
        end
    end
end