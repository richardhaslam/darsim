%%  Basis functions updater
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 16 August 2016
%Last modified: 09 August 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef bf_updater_ms_geothermal < bf_updater_ms
    properties
    end
    methods
        function ConstructPressureSystem(obj, ProductionSystem, FluidModel, FineGrid, CrossConnections)
            % Builds fine-scale incompressible pressure system.
            % In this function, mobility (Mob) is obtained based on
            % viscosity as no saturation is involved.
            K = ProductionSystem.Reservoir.K;
            Mob = FluidModel.ComputePhaseMobilities(ProductionSystem.Reservoir.State);
            obj.A = obj.MediumPressureSystem(FineGrid, K, Mob);
            obj.AddWellsToPressureMatrix(ProductionSystem.Wells, K, Mob, FineGrid.N)
        end
    end
end