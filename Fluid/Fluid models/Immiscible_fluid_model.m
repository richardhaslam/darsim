% Immiscible Fluid model base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 14 July 2016
%Last modified: 14 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Immiscible_fluid_model < fluid_model
    properties
    end
    methods
        function obj = Immiscible_fluid_model(n_phases)
            obj@fluid_model(n_phases, n_phases);
        end
        function Status = InitializeReservoir(obj, Status)
            % Composition in this case is fixed to be 1 and 0
            Status.x1(:,1) = 1;
            Status.x1(:,2) = 0;
            
            % Compute Phase Density
            for i=1:obj.NofPhases
                Status.rho(:, i) = obj.Phases(i).ComputeDensity(Status.p);
            end
            
            % Two phase, two component saturation updater    
            Status.S = Status.rho(:,2).*(Status.x1(:,2) - Status.z(:,1))./(Status.rho(:,1).*(Status.z(:,1)...
                - Status.x1(:,1)) + Status.rho(:,2).*(Status.x1(:,2) - Status.z(:,1))); 
        end
        function Inj = InitializeInjectors(obj, Inj)
            for i=1:length(Inj)
                Inj(i).z = [1 0];
                Inj(i).x1 = [1 0];
                Inj(i).s = 1;
                for ph=1:obj.NofPhases
                    Inj(i).rho(:, ph)= obj.Phases(ph).ComputeDensity(Inj(i).p);
                end
                Inj(i).x2 = 1 - Inj(i).x1;
                Inj(i).Mob = obj.ComputePhaseMobilities(Inj(i).s);   
            end
        end
    end
end