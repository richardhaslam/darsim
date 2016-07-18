% Immiscible Formulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 12 July 2016
%Last modified: 12 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Immiscible_formulation < formulation
    properties
        
    end
    methods
        function BuildResidual(obj)
        end
        function BuildJacobian(obj)
        end
        function Status = UpdateStatus(obj, delta, Status, FluidModel)
            % Update Solution
            Status.p = delta(1:end/2);
            Status.S = delta(end/2+1:end);
            
            % Remove unphysical solutions
            Status.S = min(Status.S,1);
            Status.S = max(Status.S,0);
            
            % Update density
            for i=1:FluidModel.NofPhases
                Status.rho(:, i) = FluidModel.Phases(i).ComputeDensity(Status.p);
            end
            
            % Update z
            Status.z = FluidModel.ComputeMassFractions(Status.s, Staus.x1, Status.rho);
            
            % Update total density
            Status.rhoT = FluidModel.ComputeTotalDensity(Status.s, Status.rho);
        end
        function BuildPressureSystem(obj)
        end
        function BuildTransportSystem(obj)
        end
    end
end