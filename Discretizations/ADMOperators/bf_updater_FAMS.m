%  Basis functions updater
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Mousa HosseiniMehr & Matteo Cusini
%TU Delft
%Created: 07 August 2017
%Last modified: 07 August 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef bf_updater_FAMS < bf_updater_ms
    properties
    end
    methods
        function ConstructPressureSystem(obj, ProductionSystem, FluidModel, ReservoirGrid, FracturesGrid, CrossConnections)
            % Reservoir
            Km = ProductionSystem.Reservoir.K;
            S = ProductionSystem.CreateGlobalVariables(ReservoirGrid, FracturesGrid, FluidModel.NofPhases, 'S_');
            Mob = FluidModel.ComputePhaseMobilities(S(:,1));
            Start = 1;
            End = ReservoirGrid.N;
            obj.A = MediumPressureSystem(obj, ReservoirGrid, Km, Mob(Start:End));
            % Fractures
            for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                Start = End + 1;
                End = Start + FracturesGrid.N(f) - 1;
                Kf = ProductionSystem.FracturesNetwork.Fractures(f).K;
                Af = MediumPressureSystem(obj, FracturesGrid.Grids(f), Kf, Mob(Start:End));
                obj.A = blkdiag(obj.A, Af);
            end
            % Non-Neighboring Connections
            Mobt = sum(Mob,2);
            for c = 1:length(CrossConnections)
                T_Geo = CrossConnections(c).T_Geo;
                i = c + ReservoirGrid.N;
                j = CrossConnections(c).Cells;
                obj.A(i,j) = obj.A(i,j) - T_Geo' .* Mobt(j)';
                obj.A(i,i) = obj.A(i,i) + sum(T_Geo.* Mobt(j));
                obj.A(j,i) = obj.A(j,i) - T_Geo.* Mobt(i);
                obj.A(sub2ind(size(obj.A), j, j)) = obj.A(sub2ind(size(obj.A), j, j)) + T_Geo.* Mobt(i);
            end
        end
        
    end
end