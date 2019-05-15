% Geothermal Fluid model base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Rhadityo ...
%TU Delft
%Created: 24 January 2018
%Last modified: 24 January 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Geothermal_2T_fluid_model < fluid_model
    properties
        AveragedTemperature
    end
    methods
        function obj = Geothermal_2T_fluid_model()
            obj@fluid_model(1, 1);
            obj.name = 'Geothermal_2T';
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
                    mu = obj.Phases(ph).ComputeViscosity(Inj(i).T);   
                end
                Inj(i).Mob = 1/mu;   
            end
        end
        function ComputePhaseDensities(obj, Status)
            % here you decide how to compute densities as function of P&T
            rho = Status.Properties('rho_1'); 
            rho.Value = obj.Phases(1).ComputeDensity(Status.Properties('P_1').Value, Status.Properties('Tf').Value);
        end
        function ComputePhaseEnthalpies(obj, Status)
            % here you decide how to compute enthalpy as function of P&T
            h = Status.Properties('h_1'); 
            h.Value = obj.Phases(1).ComputeEnthalpy(Status.Properties('P_1').Value, Status.Properties('Tf').Value);
        end
        function ComputePhaseViscosities(obj, Status)
            % here you decide how to compute viscosity as function of T
            mu = Status.Properties('mu_1'); 
            mu.Value = obj.Phases(1).ComputeViscosity(Status.Properties('Tf').Value);
        end
        function dmu = ComputeDmu(obj, Status)
            [~,dmu] = obj.Phases.ComputeViscosity(Status.Properties('Tf').Value);
%             dmu = obj.Phases.ComputeDmu(Status.Properties('T').Value);
        end
        function Mob = ComputePhaseMobilities(obj, mu) 
            Mob = 1./mu; 
        end
        function dMob = ComputeDMobdT(obj, Status) 
            % fill derivative of d(1/mu)/dT
            mu = Status.Properties('mu_1').Value;
            dmu = ComputeDmu(obj, Status);
            dMob = -dmu./(mu.^2);
        end
        function drho = ComputeDrho(obj, Status)
            drhodp = obj.Phases.ComputeDrhoDp(Status.Properties('P_1').Value, Status.Properties('Tf').Value);
            drhodT = obj.Phases.ComputeDrhoDT(Status.Properties('P_1').Value, Status.Properties('Tf').Value);
            drho = [drhodp, drhodT];
        end
        function dh = ComputeDh(obj, Status)
            dhdp = obj.Phases.ComputeDhDp(Status.Properties('P_1').Value, Status.Properties('Tf').Value);
            dhdT = obj.Phases.ComputeDhDT(Status.Properties('P_1').Value, Status.Properties('Tf').Value);
            dh = [dhdp, dhdT];
        end
        function v = ComputeVelocity(obj, Reservoir, mu)
%             virtual call
        end
    end
end