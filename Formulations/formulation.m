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
            obj.UpWind =  struct('x',{},'y',{},'z',{});
            obj.U =  struct('x',{},'y',{},'z',{});
        end
        function SavePhaseState(obj)
        end
        function Reset(obj)    
        end
        function UpWindAndPhaseRockFluxes(obj, Grid, Phases, Status)
            obj.GravityModel.ComputeInterfaceDensities(Grid.Nx, Grid.Ny, Grid.Nz, [Status.Properties('rho_1').Value, Status.Properties('rho_2').Value]);
            % Compute phase rock velocities and Upwind operators
            for i=1:obj.NofPhases
                P = Status.Properties(['P_', num2str(i)]).Value;
                [obj.UpWind(i), obj.U(i)] = Phases(i).UpWindAndRockFluxes(Grid, P, obj.GravityModel.RhoInt(i));
            end
        end
    end
end