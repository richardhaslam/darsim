% Immiscible Fluid model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini and Barnaby Fryer
%TU Delft
%Created: 14 July 2016
%Last modified: 28 September 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef BO_fluid_model < Comp_fluid_model
    properties
        Pref = 1e5; % Pref for Rs computation
        Rs
        dRs
    end
    methods
        function obj = BO_fluid_model(n_phases, n_comp)
            obj@Comp_fluid_model(n_phases, n_comp);
            obj.name = 'Black Oil';
        end
%         function SinglePhase = Flash(obj, Status)
%             % Define SinglePhase objects
%             SinglePhase.onlyliquid = zeros(length(Status.p), 1);
%             SinglePhase.onlyvapor = zeros(length(Status.p), 1);
%             
%             % Rs
%             obj.ComputeRs(Status);
%             
%             Status.x1(:,1) = 1;
%             Status.x1(:, 2) = obj.Components(1).rho * obj.Rs(:,2) ./ (obj.Components(2).rho * ones(length(obj.Rs), 1) + obj.Components(1).rho * obj.Rs(:, 2));
%             
%             %Recognize single phase cells and fix their xs to be equal to z
%             SinglePhase.onlyliquid(Status.x1(:, 2) >= Status.z(:,1)) = 1;
%             Status.x1(SinglePhase.onlyliquid == 1, 2) = Status.z(SinglePhase.onlyliquid == 1, 1);
%         end
        function InitializeInjectors(obj, Inj)
            % Loop over all injectors
            for i=1:length(Inj)
                Inj(i).z = [1 0];
                Inj(i).x1 = [1 0 ];                
                for ph=1:obj.NofPhases
                    Inj(i).rho(:, ph) = obj.Phases(ph).ComputeDensity(Inj(i), obj.Components, 0);
                end
                Inj(i).S = 1;
                Inj(i).x2 = 1 - Inj(i).x1;
                Inj(i).Mob = obj.ComputePhaseMobilities(Inj(i).S);   
            end            
        end
        function ComputePhaseDensities(obj, Status)
            obj.ComputeRs(Status);            
            for i=1:obj.NofPhases
                Status.rho(:, i) = obj.Phases(i).ComputeDensity(Status, obj.Components, obj.Rs(:,i));
            end
        end
        function k = ComputeKvalues(obj, Status)
            obj.Rs = zeros(length(Status.p), 2);
            obj.ComputeRs(Status);
            k = obj.KvaluesCalculator.Compute(obj.Components, obj.Rs);
        end
        function dkdp = DKvalDp(obj, p)
            dkdp = obj.KvaluesCalculator.DKvalDp(obj.Components, obj.Rs, obj.dRs);
        end
        function ComputeRs(obj, Status)
            Pdim = Status.p/obj.Pref;
            for i=1:obj.NofPhases
                [obj.Rs(:,i), obj.dRs(:,i)] = obj.Phases(i).ComputeRs(Pdim);
                obj.dRs = obj.dRs / obj.Pref;
            end
        end
        function drho = DrhoDp(obj, p)
            drho = zeros(length(p), obj.NofPhases);
            for i=1:obj.NofPhases
                drho(:, i) = obj.Phases(i).DrhoDp(p, obj.Components, obj.dRs(:,i));
            end
        end
    end
end