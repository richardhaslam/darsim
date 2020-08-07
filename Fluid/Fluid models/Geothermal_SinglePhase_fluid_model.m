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
        TablePH
        TablePT
        Ptable = (10:1:200).*1e5;  %[Pa]
        Ttable = (274:1:625);      %[K]

        Pindex % The index of pressure value for lookup in the tables
        Hindex % The index of enthalpy value for lookup in the tables
        Pstepsize = 1e5;
        Tstepsize = 1; % implement this
        
        Pgrid
        Tgrid
%         p
%         T

    end
    methods
        function obj = Geothermal_SinglePhase_fluid_model()
            obj@fluid_model(1, 1);
            obj.name = 'Geothermal_SinglePhase';
        end
        function SinglePhase = Flash(obj, Status)
            SinglePhase (:) = 1;
        end
        function InitializeInjectors(obj, Inj)
            % Assuming saturation of injection phase is 1.0. We are assuming we are only injecting 1 phase, i.e. water
            for i=1:length(Inj)
                Inj(i).z = 1;
%                 Inj(i).x = [1 0];
%                 Inj(i).S = 1;
                Inj(i).Mob = zeros(length(Inj(i).p),obj.NofPhases);
                
                if strcmp(Inj(i).BC_Formulation, 'Temperature')
                    Inj(i).h(:, 1) = obj.Phases(1).ComputeWaterEnthalpy(Inj(i).p, Inj(i).T);
                    Inj(i).rho(:, 1)= obj.Phases(1).ComputeWaterDensity(Inj(i).p, Inj(i).T);
                    mu = obj.Phases(1).ComputeWaterViscosity(Inj(i).T);
                elseif strcmp(Inj(i).BC_Formulation, 'Enthalpy')
%                     Inj(i).S = obj.Phases(1).GetSaturation(obj.Pgrid, obj.Hgrid, obj.TablePH.S_1, Inj(i).h(:,1), Inj(i).p);
                    Inj(i).T = interp2(obj.Hgrid, obj.Pgrid, obj.TablePH.Temperature, Inj(i).h(:,1), Inj(i).p, 'linear');
                    Inj(i).rho(:, 1) = obj.Phases(1).GetDensity(obj.Pgrid, obj.Hgrid, obj.TablePH.rho_1, Inj(i).h(:,1), Inj(i).p);
                    mu = obj.Phases(1).GetViscosity(obj.Pgrid, obj.Hgrid, obj.TablePH.mu_1, Inj(i).h(:,1), Inj(i).p);
                end
                    
                Inj(i).h(:, 2:obj.NofPhases) = 0;
                Inj(i).rho(:, 2:obj.NofPhases) = 0; 
                Inj(i).Mob(:, 1) = 1/mu; % injecting only water means rel.perm of water is 1.0 
                Inj(i).Mob(:, 2:obj.NofPhases) = 0;            
                
%                 Inj(i).h(:, 2:obj.NofPhases) = Inj(i).h(:,1);
%                 Inj(i).rho(:, 2:obj.NofPhases) = Inj(i).rho(: ,1); 
%                 Inj(i).Mob(:, 1) = 1/mu; % injecting only water means rel.perm of water is 1.0 
%                 Inj(i).Mob(:, 2:obj.NofPhases) = Inj(i).Mob(:, 1);            
            end
        end
        
%         function ComputePhaseDensities(obj, Status)
%             % here you decide how to compute densities as function of P&T
%             rho = Status.Properties('rho_1'); 
%             rho.Value = obj.Phases(1).ComputeDensity(Status.Properties('P_1').Value, Status.Properties('T').Value);
%         end
        function AddPhaseConductivities(obj, Status)
            cond = Status.Properties('cond_1');
            cond.Value = obj.Phases(1).AddConductivity(Status.Properties('P_1').Value, Status.Properties('T').Value);
        end
%         function ComputePhaseEnthalpies(obj, Status)
%             % here you decide how to compute enthalpy as function of P&T
%             h = Status.Properties('h_1'); 
%             h.Value = obj.Phases(1).ComputeEnthalpy(Status.Properties('P_1').Value, Status.Properties('T').Value);
%         end
%         function ComputePhaseViscosities(obj, Status)
%             mu = Status.Properties('mu_1'); 
%             [mu.Value,~,~] = obj.Phases(1).ComputeViscosity(Status.Properties('T').Value);
%         end
%         function [dmudT,d2mudT2] = ComputeDmuDT(obj, Status)
%                 [~,dmudT,d2mudT2] = obj.Phases.ComputeViscosity(Status.Properties('T').Value);
%         end
%         function Mob = ComputePhaseMobilities(obj, mu) 
%             Mob = 1./mu;
%         end
%         function dMobdp = ComputeDMobDp(obj,status)
%             % For now, mobility has no dependency on pressure.
%         end
%         function [dMobdT,d2MobdT2] = ComputeDMobdT(obj, Status) 
%             mu = Status.Properties('mu_1').Value;
%             [dmudT,d2mudT2] = obj.ComputeDmuDT(Status);
%             dMobdT = -dmudT./(mu.^2);
%             d2MobdT2 = (2.*dmudT.^2./mu.^3) - d2mudT2./mu.^2;
%         end
%         function drhodp = ComputeDrhoDp(obj, Status)
%             drhodp = obj.Phases.ComputeDrhoDp(Status.Properties('P_1').Value, Status.Properties('T').Value);
%         end
%         function [drhodT,d2rhodT2] = ComputeDrhoDT(obj, Status)
%             [drhodT,d2rhodT2] = obj.Phases.ComputeDrhoDT(Status.Properties('P_1').Value, Status.Properties('T').Value);
%         end
%         function dhdp = ComputeDhDp(obj, Status)
%             dhdp = obj.Phases.ComputeDhDp(Status.Properties('P_1').Value, Status.Properties('T').Value);
%         end
%         function [dhdT,d2hdT2] = ComputeDhDT(obj, Status)
%             [dhdT,d2hdT2] = obj.Phases.ComputeDhDT(Status.Properties('P_1').Value, Status.Properties('T').Value);
%         end


        function ComputeTableGrid(obj)
            [obj.Tgrid,obj.Pgrid] = meshgrid(obj.Ttable,obj.Ptable);
        end
%         function GetPTValues(obj, Status)
%             obj.p = Status.Properties('P_1').Value;
%             obj.T = Status.Properties('T').Value;
%         end
      
        % Linear interpolation
        function ComputePhaseDensities(obj, Status)
            p = Status.Properties('P_1').Value;
            T = Status.Properties('T').Value;
            for i=1:obj.NofPhases
                rho = Status.Properties(['rho_', num2str(i)]);
                rhoTable = obj.TablePT.('rho');
                rho.Value = obj.Phases(i).ComputeDensity(obj.Pgrid, obj.Tgrid, rhoTable, T, p);
            end
        end
%         function AddPhaseConductivities(obj, Status) 
%             p = Status.Properties('P_1').Value;
%             T = Status.Properties('T').Value;
%             for i=1:obj.NofPhases
%                 ThermCond = Status.Properties(['cond_',num2str(i)]);
%                 ThermCondTable = obj.TablePT.('cond');
%                 ThermCond.Value = obj.Phases(i).AddConductivity(obj.Pgrid, obj.Tgrid, ThermCondTable, T, p);
%             end
%         end
        function ComputePhaseViscosities(obj, Status)
            p = Status.Properties('P_1').Value;
            T = Status.Properties('T').Value;
            for i=1:obj.NofPhases
                mu = Status.Properties(['mu_',num2str(i)]);
                muTable = obj.TablePT.('mu');
                mu.Value = obj.Phases(i).ComputeViscosity(obj.Pgrid, obj.Tgrid, muTable, T, p);
            end
        end
        function ComputePhaseEnthalpies(obj, Status)
            % here you decide how to compute enthalpy as function of P&T
            p = Status.Properties('P_1').Value;
            T = Status.Properties('T').Value;
            H = Status.Properties('h_1'); 
            PhaseEnthalpyTable = obj.TablePT.('h');
            H.Value = obj.Phases(1).ComputeEnthalpy(obj.Pgrid, obj.Tgrid, PhaseEnthalpyTable, T, p);
        end

        %%%
        function [dmudT,d2mudT2] = ComputeDmuDT(obj, Status)
            p = Status.Properties('P_1').Value;
            T = Status.Properties('T').Value;
            muTable = obj.TablePT.('mu');
            [dmudT,d2mudT2] = obj.Phases(1).ComputeDmuDT(obj.Pgrid, obj.Tgrid, muTable, T, p);              
        end
        function Mob = ComputePhaseMobilities(obj, mu) 
            Mob = 1./mu;
        end
        function dMobdp = ComputeDMobDp(obj,status)
            % For now, mobility has no dependency on pressure.
        end
        function [dMobdT,d2MobdT2] = ComputeDMobdT(obj, Status)
            p = Status.Properties('P_1').Value;
            T = Status.Properties('T').Value;
            mu = Status.Properties('mu_1').Value;
            muTable = obj.TablePT.('mu');
            [dmudT,d2mudT2] = obj.Phases(1).ComputeDmuDT(obj.Pgrid, obj.Tgrid, muTable, T, p);              
            dMobdT = -dmudT./(mu.^2);
            d2MobdT2 = (2.*dmudT.^2./mu.^3) - d2mudT2./mu.^2;
        end
        function drhodp = ComputeDrhoDp(obj, Status)
            p = Status.Properties('P_1').Value;
            T = Status.Properties('T').Value;
            rhoTable = obj.TablePT.('rho');
            drhodp = obj.Phases(1).ComputeDrhoDp(obj.Pgrid, obj.Tgrid, rhoTable, T, p);
        end
        function [drhodT,d2rhodT2] = ComputeDrhoDT(obj, Status)
            p = Status.Properties('P_1').Value;
            T = Status.Properties('T').Value;
            rhoTable = obj.TablePT.('rho');
            [drhodT,d2rhodT2] = obj.Phases.ComputeDrhoDT(obj.Pgrid, obj.Tgrid, rhoTable, T, p);
        end       
        function dhdp = ComputeDhDp(obj, Status)
            p = Status.Properties('P_1').Value;
            T = Status.Properties('T').Value;
            PhaseEnthalpyTable = obj.TablePT.('h');
            dhdp = obj.Phases.ComputeDhDp(obj.Pgrid, obj.Tgrid, PhaseEnthalpyTable, T, p);
        end       
        function [dhdT,d2hdT2] = ComputeDhDT(obj, Status)
            p = Status.Properties('P_1').Value;
            T = Status.Properties('T').Value;
            PhaseEnthalpyTable = obj.TablePT.('h');
            [dhdT,d2hdT2] = obj.Phases.ComputeDhDT(obj.Pgrid, obj.Tgrid, PhaseEnthalpyTable, T, p);
        end
        
        
        function ComputeThermalConductivity(obj, Reservoir)
            % (1-phi)*C_r + phi*Sw*C_w + phi*Ss*C_s : vector
            CondEff = Reservoir.State.Properties('CondEff');
            CondEff.Value = (1 - Reservoir.Por) .* Reservoir.K_Cond_rock;
            for i=1:obj.NofPhases
                cond = Reservoir.State.Properties(['cond_',num2str(i)]);
                S = 1;
                CondEff.Value = CondEff.Value + Reservoir.Por .* cond.Value .* S;
            end
        end
        function v = ComputeVelocity(obj, Reservoir, mu)
%             virtual call
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         function obj = Geothermal_SinglePhase_fluid_model()
%             obj@fluid_model(1, 1);
%             obj.name = 'Geothermal_SinglePhase';
%         end
%         function SinglePhase = Flash(obj, Status)
%             SinglePhase (:) = 1;
%         end
%         function InitializeInjectors(obj, Inj)
%             for i=1:length(Inj)
%                 Inj(i).z = 1;
%                 Inj(i).x = [1 0];
%                 Inj(i).S = 1;
%                 for ph=1:obj.NofPhases
%                     Inj(i).rho(:, ph)= obj.Phases(ph).ComputeDensity(Inj(i).p, Inj(i).T);
%                     Inj(i).h(:, ph)= obj.Phases(ph).ComputeEnthalpy(Inj(i).p, Inj(i).T);
%                     mu = obj.Phases(ph).ComputeViscosity(Inj(i).T);   
%                 end
%                 Inj(i).Mob = 1/mu;   
%             end
%         end
%         function ComputePhaseDensities(obj, Status)
%             % here you decide how to compute densities as function of P&T
%             rho = Status.Properties('rho_1'); 
%             rho.Value = obj.Phases(1).ComputeDensity(Status.Properties('P_1').Value, Status.Properties('T').Value);
%         end
%         function AddPhaseConductivities(obj, Status)
%             cond = Status.Properties('cond_1');
%             cond.Value = obj.Phases(1).AddConductivity(Status.Properties('P_1').Value, Status.Properties('T').Value);
%         end
%         function ComputePhaseEnthalpies(obj, Status)
%             % here you decide how to compute enthalpy as function of P&T
%             h = Status.Properties('h_1'); 
%             h.Value = obj.Phases(1).ComputeEnthalpy(Status.Properties('P_1').Value, Status.Properties('T').Value);
%         end
%         function ComputePhaseViscosities(obj, Status)
%             mu = Status.Properties('mu_1'); 
%             [mu.Value,~,~] = obj.Phases(1).ComputeViscosity(Status.Properties('T').Value);
%         end
%         function [dmudT,d2mudT2] = ComputeDmuDT(obj, Status)
%                 [~,dmudT,d2mudT2] = obj.Phases.ComputeViscosity(Status.Properties('T').Value);
%         end
%         function Mob = ComputePhaseMobilities(obj, mu) 
%             Mob = 1./mu;
%         end
%         function dMobdp = ComputeDMobDp(obj,status)
%             % For now, mobility has no dependency on pressure.
%         end
%         function [dMobdT,d2MobdT2] = ComputeDMobdT(obj, Status) 
%             mu = Status.Properties('mu_1').Value;
%             [dmudT,d2mudT2] = obj.ComputeDmuDT(Status);
%             dMobdT = -dmudT./(mu.^2);
%             d2MobdT2 = (2.*dmudT.^2./mu.^3) - d2mudT2./mu.^2;
%         end
%         function drhodp = ComputeDrhoDp(obj, Status)
%             drhodp = obj.Phases.ComputeDrhoDp(Status.Properties('P_1').Value, Status.Properties('T').Value);
%         end
%         function [drhodT,d2rhodT2] = ComputeDrhoDT(obj, Status)
%             [drhodT,d2rhodT2] = obj.Phases.ComputeDrhoDT(Status.Properties('P_1').Value, Status.Properties('T').Value);
%         end
%         function dhdp = ComputeDhDp(obj, Status)
%             dhdp = obj.Phases.ComputeDhDp(Status.Properties('P_1').Value, Status.Properties('T').Value);
%         end
%         function [dhdT,d2hdT2] = ComputeDhDT(obj, Status)
%             [dhdT,d2hdT2] = obj.Phases.ComputeDhDT(Status.Properties('P_1').Value, Status.Properties('T').Value);
%         end
%         function v = ComputeVelocity(obj, Reservoir, mu)
% %             virtual call
%         end

    end
end