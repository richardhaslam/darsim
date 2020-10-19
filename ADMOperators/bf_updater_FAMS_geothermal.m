%%  Basis functions updater F-AMS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Mousa HosseiniMehr & Matteo Cusini
%TU Delft
%Created: 07 August 2017
%Last modified: 14 August 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef bf_updater_FAMS_geothermal < bf_updater_FAMS
    properties
    end
    methods
        function ConstructPressureSystem(obj, ProductionSystem, FluidModel, FineGrid, CrossConnections)
            % Builds fine-scale incompressible pressure system.
            % In this function, mobility (Mob) is obtained based on
            % viscosity as no saturation is involved.
            % Reservoir
            Km = ProductionSystem.Reservoir.K;
            ReservoirMobility = FluidModel.ComputePhaseMobilities(ProductionSystem.Reservoir.State);
            obj.Amedia{1} = obj.MediumPressureSystem(FineGrid(1), Km, ReservoirMobility);
            obj.A = obj.Amedia{1};
            obj.AddWellsToPressureMatrix(ProductionSystem.Wells, Km, ReservoirMobility, FineGrid(1).N);
            GlobalMobility = ReservoirMobility;
            % Fractures
            for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                Kf = ProductionSystem.FracturesNetwork.Fractures(f).K;
                FractureMobility = FluidModel.ComputePhaseMobilities(ProductionSystem.FracturesNetwork.Fractures(f).State);
                obj.Amedia{1+f} = obj.MediumPressureSystem(FineGrid(1+f), Kf, FractureMobility);
                obj.A = blkdiag(obj.A, obj.Amedia{1+f});
                GlobalMobility = [GlobalMobility; FractureMobility];
            end
            % Non-Neighboring Connections
            TotalMobility = sum(GlobalMobility,2);
            for c = 1:length(CrossConnections)
                T_Geo = CrossConnections(c).T_Geo;
                i = c + FineGrid(1).N;
                j = CrossConnections(c).Cells;
                obj.A(i,j) = obj.A(i,j) - T_Geo' .* TotalMobility(j)';
                obj.A(i,i) = obj.A(i,i) + sum(T_Geo.* TotalMobility(j));
                obj.A(j,i) = obj.A(j,i) - T_Geo.* TotalMobility(i);
                obj.A(sub2ind(size(obj.A), j, j)) = obj.A(sub2ind(size(obj.A), j, j)) + T_Geo.* TotalMobility(i);
            end 
        end
    end
end