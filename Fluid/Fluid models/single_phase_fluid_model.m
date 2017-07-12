% SinglePhase Fluid model base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 2 March 2017
%Last modified: 4 March 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef single_phase_fluid_model < fluid_model
    properties
        
    end
    methods
        function obj = single_phase_fluid_model()
            obj@fluid_model(1, 1);
            obj.name = 'SinglePhase';
        end
        function SinglePhase = Flash(obj, Status)
            SinglePhase (:) = 1;
        end
        function InitializeInjectors(obj, Inj)
            for i=1:length(Inj)
                Inj(i).z = 1;
                Inj(i).x = [1 0];
                Inj(i).S = 1;
                for ph=1:obj.NofPhases
                    Inj(i).rho(:, ph)= obj.Phases(ph).ComputeDensity(Inj(i).p);
                end
                Inj(i).Mob = 1/obj.Phases(1).mu;   
            end
        end
        function ComputePhaseDensities(obj, Status)
            rho = Status.Properties('rho_1'); 
            rho.Value = obj.Phases(1).ComputeDensity(Status.Properties('P_1').Value, obj.Components);
        end
         function ComputeTotalDensity(obj, Status)
            % Compute the total density
            rhoT = Status.Properties('rhoT');
            rhoT.Value = Status.Properties('rho_1').Value; 
            % For 1 phase rhoT is rho1
         end
        function Mob = ComputePhaseMobilities(obj, s)
            Mob = ones(length(s), obj.NofPhases);
            Mob(:,1) = 1/obj.Phases(1).mu;
        end
         function dMob = DMobDS(obj, S)
            dMob = zeros(length(S), 1);
        end
        function dPc = DPcDS(obj, S)
            dPc = zeros(length(S), 1);
        end
    end
end