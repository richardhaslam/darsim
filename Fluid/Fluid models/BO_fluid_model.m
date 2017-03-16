% Immiscible Fluid model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini and Barnaby Fryer
%TU Delft
%Created: 14 July 2016
%Last modified: 7 March 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef BO_fluid_model < Comp_fluid_model
    properties
        Rs
        dRs
    end
    methods
        function obj = BO_fluid_model(n_phases, n_comp)
            obj@Comp_fluid_model(n_phases, n_comp);
            obj.name = 'Black Oil';
        end
        function SinglePhase = Flash(obj, Status)
           SinglePhase = obj.FlashCalculator.Flash(Status, obj.Components, obj.Phases);
        end
        function InitializeInjectors(obj, Inj)
            % Loop over all injectors
            for i=1:length(Inj)
                Inj(i).z = [1 0];
                Inj(i).x = [1 0 0 1];        
                for ph=1:obj.NofPhases
                    Inj(i).rho(:, ph) = obj.Phases(ph).ComputeDensity(Inj(i).p, obj.Components, zeros(length(Inj(i).Cells), 1));
                end
                Inj(i).S = 1;
                Inj(i).Mob = obj.ComputePhaseMobilities(Inj(i).S);   
            end            
        end
        function ComputePhaseDensities(obj, Status, SinglePhase)
            [obj.Rs, obj.dRs] = obj.FlashCalculator.KvaluesCalculator.ComputeRs(Status, obj.Phases);
            for i=1:obj.NofPhases
                % Correct Rs and dRs for undersaturated cells
                [obj.Rs(:,i), obj.dRs(:,i)] = obj.Phases(i).RsOfUnderSaturatedPhase(Status.Properties('z_1').Value, obj.Components, obj.Rs(:,i), obj.dRs(:,i), SinglePhase);
                rho = Status.Properties(['rho_', num2str(i)]);
                rho.Value = obj.Phases(i).ComputeDensity(Status.Properties('P_2').Value, obj.Components, obj.Rs(:,i));
            end
        end
        function k = ComputeKvalues(obj, Status)            
            k = obj.FlashCalculator.KvaluesCalculator.Compute(Status, obj.Components, obj.Phases);
        end
        function dkdp = DKvalDp(obj, Status)
            dkdp = obj.FlashCalculator.KvaluesCalculator.DKvalDp(Status, obj.Components, obj.Phases);
        end
        function drho = DrhoDp(obj, p, SinglePhase)
            drho = zeros(length(p), obj.NofPhases);
            for i=1:obj.NofPhases
                drho(:, i) = obj.Phases(i).DrhoDp(p, obj.Components, obj.Rs(:,i), obj.dRs(:,i));
            end
        end
        function drho = DrhoDz(obj, Status, SinglePhase)
            N = length(Status.Properties('P_2').Value);
            p = Status.Properties('P_2').Value;
            z = Status.Properties('z_1').Value;
            drho = zeros(N, obj.NofPhases);
            drho(:, 2) = obj.Phases(2).DrhoDz(p, z, obj.Components, SinglePhase);
        end
    end
end