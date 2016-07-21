% Immiscible Fluid model base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 14 July 2016
%Last modified: 18 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Immiscible_fluid_model < fluid_model
    properties
    end
    methods
        function obj = Immiscible_fluid_model(n_phases)
            obj@fluid_model(n_phases, n_phases);
        end
        function InitializeReservoir(obj, Status)
            % Define initial values
            P_init = linspace(1e5, 10e4, length(Status.p));
            z_init = ones(length(Status.p), 1)*0.0;
            
            % Assign initial valus
            Status.p = Status.p .* P_init';
            Status.z(:,1) = z_init;
            Status.z(:,2) = 1 - z_init;
            
            % Composition in this case is fixed to be 1 and 0
            Status.x1(:,1) = 1;
            Status.x1(:,2) = 0;
            
            % Compute Phase Density
            for i=1:obj.NofPhases
                Status.rho(:, i) = obj.Phases(i).ComputeDensity(Status.p);
            end
            
            SinglePhase.onlyvapor (Status.z(:,1) == 1) = 1;
            SinglePhase.onlyliquid (Status.z(:,2) == 1) = 1;
            
            % Two phase, two component saturation updater    
            obj.ComputePhaseSaturation(Status, SinglePhase);
        end
        function InitializeInjectors(obj, Inj)
            for i=1:length(Inj)
                Inj(i).z = [1 0];
                Inj(i).x1 = [1 0];
                Inj(i).S = 1;
                for ph=1:obj.NofPhases
                    Inj(i).rho(:, ph)= obj.Phases(ph).ComputeDensity(Inj(i).p);
                end
                Inj(i).x2 = 1 - Inj(i).x1;
                Inj(i).Mob = obj.ComputePhaseMobilities(Inj(i).S);   
            end
        end
    end
end