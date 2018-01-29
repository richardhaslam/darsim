% Geothermal Fluid model base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Rhadityo ...
%TU Delft
%Created: 24 January 2018
%Last modified: 24 January 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Geothermal_fluid_model < fluid_model
    properties
        
    end
    methods
        function obj = Geothermal_fluid_model()
            obj@fluid_model(1, 1);
            obj.name = 'Geothermal';
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
                    Inj(i).rho(:, ph)= obj.Phases(ph).ComputeDensity(Inj(i).p, Inj(i).T);
                    Inj(i).h(:, ph)= obj.Phases(ph).ComputeEnthalpy(Inj(i).p, Inj(i).T);
                    mu = obj.Phases(ph).ComputeViscosity(Inj(i).p, Inj(i).T);   
                end
                Inj(i).Mob = 1/mu;   
            end
        end
        function ComputePhaseDensities(obj, Status)
            % here you decide how to compute densities
            rho = Status.Properties('rho_1'); 
            rho.Value = obj.Phases(1).ComputeDensity(Status.Properties('P_1').Value, Status.Properties('T').Value);
        end
        function ComputePhaseEnthalpies(obj, Status)
            % here you decide how to compute densities
            h = Status.Properties('h_1'); 
            h.Value = obj.Phases(1).ComputeEnthalpy(Status.Properties('P_1').Value, Status.Properties('T').Value);
        end
        function ComputePhaseViscosities(obj, Status)
            % here you decide how to compute densities
            mu = Status.Properties('mu_1'); 
            mu.Value = obj.Phases(1).ComputeDensity(Status.Properties('P_1').Value, Status.Properties('T').Value);
        end
        function Mob = ComputePhaseMobilities(obj, Status)
            % call computation of viscosity first
            mu = Status.Properties('mu_1'); 
            Mob = 1/mu.Value;
        end
        function dMob = DMobDp(obj, Status)
            % fill it in
        end
        function dMob = DMobDT(obj, Status)
           % fill it in
        end
    end
end