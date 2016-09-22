%  FIM Formulation base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 12 July 2016
%Last modified: 19 September 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef fim_formulation < handle
    properties
        UpWindPh1
        UpWindPh2
        Uph1
        Uph2
        Tph1
        Tph2
        Gph1
        Gph2
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
        function UpWindAndPhaseRockFluxes(obj, Grid, Phases, Status)
            
            obj.GravityModel.ComputeInterfaceDensities(Grid.Nx, Grid.Ny, Status.rho);
            % Compute phase rock velocities and Upwind operators
            [obj.UpWindPh1, obj.Uph1] = Phases(1).UpWindAndRockFluxes(Grid, Status.p-Status.pc, obj.GravityModel.RhoInt1, Status.rho(:,1));
            [obj.UpWindPh2, obj.Uph2] = Phases(2).UpWindAndRockFluxes(Grid, Status.p, obj.GravityModel.RhoInt2, Status.rho(:,2));
        end
    end
end