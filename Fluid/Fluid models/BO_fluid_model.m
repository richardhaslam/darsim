% Immiscible Fluid model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 14 July 2016
%Last modified: 28 September 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef BO_fluid_model < Comp_fluid_model
    properties
        Pref % Pref for Rs computation
        Rs
        dRs
    end
    methods
        function obj = BO_fluid_model(n_phases, n_comp)
            obj@Comp_fluid_model(n_phases, n_comp);
            obj.Pref = 1e7;
        end
        function SinglePhase = Flash(obj, Status)
            % Define SinglePhase objects
            SinglePhase.onlyliquid = zeros(length(Status.p), 1);
            SinglePhase.onlyvapor = zeros(length(Status.p), 1);
            
            %% - Dimensionless pressure!
            Pdim = Status.p / obj.Pref;
            
            k = obj.ComputeKvalues(Pdim, Status.T);
           
            % Define compositions... 
            
            %Recognize single phase cells and fix their xs to be equal to z
            SinglePhase.onlyliquid(Status.x1(:, 2) >= Status.z(:,1)) = 1;
            Status.x1(SinglePhase.onlyliquid == 1, 2) = Status.z(SinglePhase.onlyliquid == 1, 1);
        end
        function ComputePhaseDensities(obj, Status)
            obj.ComputeRs(Status);            
            for i=1:obj.NofPhases
                Status.rho(:, i) = obj.Phases(i).ComputeDensity(Status.p, obj.Rs(i));
            end
        end
        function k = ComputeKvalues(obj, p, T)
            k = obj.KvaluesCalculator.Compute(p, T, obj.Components, obj.Rs);
        end
        function dkdp = DKvalDp(obj, p)
            dkdp = obj.DKvalDp(p, obj.Rs, obj.dRs);
        end
        function ComputeRs(obj, Status)
            Pdim = Status.p/obj.Pref;
            for i=1:objNofPhases
                [obj.Rs(:,i), obj.dRs(:,i)] = obj.Phases(i).ComputeRs();
            end
        end
    end
end