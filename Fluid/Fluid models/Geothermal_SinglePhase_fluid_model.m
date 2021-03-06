% Geothermal Fluid model base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Rhadityo ...
%TU Delft
%Created: 24 January 2018
%Last modified: 24 January 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Geothermal_SinglePhase_fluid_model < fluid_model
    properties
    end
    methods
        function obj = Geothermal_SinglePhase_fluid_model()
            obj@fluid_model(1, 1);
            obj.name = 'Geothermal_SinglePhase';
        end
        function SinglePhase = Flash(obj, Status)
            SinglePhase (:) = 1;
        end
        function ComputePhaseDensities(obj, Status)
            % here you decide how to compute densities as function of P&T
            rho = Status.Properties('rho_1'); 
            rho.Value = obj.Phases(1).ComputeDensityBasedOnTemperature(Status.Properties('P_1').Value, Status.Properties('T').Value);
        end
        function AddPhaseConductivities(obj, Status)
            cond = Status.Properties('cond_1');
            cond.Value = obj.Phases(1).AddConductivity();
        end
        function ComputePhaseEnthalpies(obj, Status)
            % Here you decide how to compute enthalpy as function of P&T
            h = Status.Properties('h_1'); 
            h.Value = obj.Phases(1).ComputeEnthalpy(Status.Properties('P_1').Value, Status.Properties('T').Value);
            hTfluid = Status.Properties('hTfluid');
            hTfluid.Value = h.Value;
        end
        function ComputePhaseViscosities(obj, Status)
            mu = Status.Properties('mu_1'); 
            [mu.Value,~,~] = obj.Phases(1).ComputeViscosity(Status.Properties('T').Value);
        end
        function [dmudT,d2mudT2] = ComputeDmuDT(obj, Status)
                [~,dmudT,d2mudT2] = obj.Phases.ComputeViscosity(Status.Properties('T').Value);
        end
        function Mob = ComputePhaseMobilities(obj, Status) 
            mu = Status.Properties('mu_1').Value;
            Mob = 1./mu;
        end
        function dMobdp = ComputeDMobDp(obj,status)
            % For now, mobility has no dependency on pressure.
        end
        function [dMobdT,d2MobdT2] = ComputeDMobdT(obj, Status) 
            mu = Status.Properties('mu_1').Value;
            [dmudT,d2mudT2] = obj.ComputeDmuDT(Status);
            dMobdT = -dmudT./(mu.^2);
            d2MobdT2 = (2.*dmudT.^2./mu.^3) - d2mudT2./mu.^2;
        end
        function ComputeThermalConductivity(obj, Medium)
            % (1-phi)*C_r + phi*Sw*C_w + phi*Ss*C_s : vector
            CondEff = Medium.State.Properties('CondEff');
            CondEff.Value = (1 - Medium.Por) .* Medium.K_Cond_rock;
            for i=1:obj.NofPhases
                cond = obj.Phases(i).AddConductivity();
                S = Medium.State.Properties(['S_',num2str(i)]);
                CondEff.Value = CondEff.Value + Medium.Por .* cond .* S.Value;
            end
        end
        function drhodp = ComputeDrhoDp(obj, Status)
            drhodp = obj.Phases.ComputeDrhoDp(Status.Properties('P_1').Value, Status.Properties('T').Value);
        end
        function [drhodT,d2rhodT2] = ComputeDrhoDT(obj, Status)
            [drhodT,d2rhodT2] = obj.Phases.ComputeDrhoDT(Status.Properties('P_1').Value, Status.Properties('T').Value);
        end
        function dhdp = ComputeDhDp(obj, Status)
            dhdp = obj.Phases.ComputeDhDp(Status.Properties('P_1').Value, Status.Properties('T').Value);
        end
        function [dhdT,d2hdT2] = ComputeDhDT(obj, Status)
            [dhdT,d2hdT2] = obj.Phases.ComputeDhDT(Status.Properties('P_1').Value, Status.Properties('T').Value);
        end
        function InitializeInjectors(obj, Inj)
            for i=1:length(Inj)
                Inj(i).z = 1;
                Inj(i).x = [1 0];
                Inj(i).S = 1;
                
                if strcmp(Inj(i).BC_Formulation, 'Temperature')
                    Inj(i).h(:, 1) = obj.Phases(1).ComputeWaterEnthalpy(Inj(i).p, Inj(i).T);
                    Inj(i).rho(:, 1)= obj.Phases(1).ComputeWaterDensity(Inj(i).p, Inj(i).T);
                    mu = obj.Phases(1).ComputeWaterViscosity(Inj(i).T);
                elseif strcmp(Inj(i).BC_Formulation, 'Enthalpy')
                    PhaseIndex = 1; % the phase index "1" refers to water phase
                    Inj(i).T = obj.Phases(1).ComputeWaterTemperature(Inj(i).p, Inj(i).h(:,1));
                    Inj(i).rho(:, 1) = obj.Phases(1).ComputeDensityBasedOnEnthalpy(1, PhaseIndex, Inj(i).p, Inj(i).h(:,1)); % the "1" is for water phase
                    mu = obj.Phases(1).ComputeWaterViscosity(Inj(i).T);
                end
                
                Inj(i).Mob = 1/mu;
            end
        end
        function v = ComputeVelocity(obj, Reservoir, mu)
%             virtual call
        end
    end
end