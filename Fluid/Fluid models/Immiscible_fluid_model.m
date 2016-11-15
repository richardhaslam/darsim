% Immiscible Fluid model base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 14 July 2016
%Last modified: 28 September 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Immiscible_fluid_model < fluid_model
    properties
        
    end
    methods
        function obj = Immiscible_fluid_model(n_phases)
            obj@fluid_model(n_phases, n_phases);
            obj.name = 'Immiscible';
        end
        function SinglePhase = Flash(obj, Status)
            % Composition in this case is fixed to be 1 and 0
            Status.x(:,1) = 1;
            Status.x(:,2) = 0;
            Status.x(:,3) = 0;
            Status.x(:,4) = 1;
            
            SinglePhase (Status.z(:,1) == 1) = 1;
            SinglePhase (Status.z(:,2) == 1) = 2;
            
        end
        function InitializeInjectors(obj, Inj)
            for i=1:length(Inj)
                Inj(i).z = [1 0];
                Inj(i).x = [1 0 0 1];
                Inj(i).S = 1;
                for ph=1:obj.NofPhases
                    Inj(i).rho(:, ph)= obj.Phases(ph).ComputeDensity(Inj(i));
                end
                Inj(i).x2 = 1 - Inj(i).x;
                Inj(i).Mob = obj.ComputePhaseMobilities(Inj(i).S);   
            end
        end
        function ComputePhaseDensities(obj, Status)
            for i=1:obj.NofPhases
                Status.rho(:, i) = obj.Phases(i).ComputeDensity(Status, obj.Components);
            end
        end
    end
end