%  FIM Formulation base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 12 July 2016
%Last modified: 26 September 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef fim_formulation < handle
    properties
        UpWind
        U
        Tph
        Gph
        Mob
        dMob
        dPc
        drho
        GravityModel
    end
    methods (Abstract)
        obj = BuildResidual(obj)
        obj = BuildJacobian(obj)
        obj = UpdateState(obj)
        obj = Reset(obj)
    end
    methods
        function obj = fim_formulation()
            obj.UpWind =  struct('x',{},'y',{});
            obj.U =  struct('x',{},'y',{});
        end
        function UpWindAndPhaseRockFluxes(obj, Grid, Phases, Status)            
            obj.GravityModel.ComputeInterfaceDensities(Grid.Nx, Grid.Ny, Status.rho);
            % Compute phase rock velocities and Upwind operators
            P(:, 2) = Status.p;
            P(:, 1) = Status.p - Status.pc;
            for i=1:2
                [obj.UpWind(i), obj.U(i)] = Phases(i).UpWindAndRockFluxes(Grid, P(:,i), obj.GravityModel.RhoInt(i));
            end
        end
    end
end