% Formulation base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 12 July 2016
%Last modified: 26 September 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef formulation < handle
    properties
        NofPhases
        UpWind
        U
        Tph
        Gph
        Mob
        dMob
        dPc
        drhodp
        GravityModel
        SinglePhase
    end
    methods (Abstract)
        obj = BuildResidual(obj)
        obj = BuildJacobian(obj)
        obj = UpdateState(obj)
    end
    methods
        function obj = formulation()
            obj.UpWind = cell(1);  %struct('x',{},'y',{},'z',{});
            obj.U =  cell(1);%struct('x',{},'y',{},'z',{});
        end
        function SavePhaseState(obj)
            % virtual call
        end
        function Reset(obj)
            % virtual call
        end
        function UpWindAndPhaseRockFluxes(obj, DiscretizationModel, Phases, ProductionSystem)
            %% Reservoir upwind
            Grid = DiscretizationModel.ReservoirGrid;
            obj.GravityModel.ComputeInterfaceDensities(Grid, ProductionSystem.Reservoir.State, 0);
            % Compute phase rock velocities and Upwind operators
            for i=1:obj.NofPhases
                P = ProductionSystem.Reservoir.State.Properties(['P_', num2str(i)]).Value;
                [obj.UpWind{i, 1}, obj.U{i, 1}] = Phases(i).UpWindAndRockFluxes(Grid, P, obj.GravityModel.RhoInt{i, 1});
            end
            %% Fractures upwind operators
            if ProductionSystem.FracturesNetwork.Active
                Grids = DiscretizationModel.FracturesGrid.Grids;
                for f=1:ProductionSystem.FracturesNetwork.NumOfFrac
                    obj.GravityModel.ComputeInterfaceDensities(Grids(f), ProductionSystem.FracturesNetwork.Fractures(f).State, f);
                    % Compute phase rock velocities and Upwind operators
                    for i=1:obj.NofPhases
                        P = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(['P_', num2str(i)]).Value;
                        [obj.UpWind{i, f+1}, obj.U{i, f+1}] = Phases(i).UpWindAndRockFluxes(Grids(f), P, obj.GravityModel.RhoInt{i, f+1});
                    end
                end
            end
        end
    end
end