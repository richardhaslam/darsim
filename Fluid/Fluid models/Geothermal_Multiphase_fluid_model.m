% Geothermal Fluid model base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Arjan Marelis
%TU Delft
%Created: 21 January 2020
%Last modified: 21 January 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Geothermal_Multiphase_fluid_model < fluid_model
    properties
        TablePH

        dSdh
    end
    methods
        function obj = Geothermal_Multiphase_fluid_model(n_phases)
            obj@fluid_model(n_phases, n_phases);
            obj.name = 'Geothermal_MultiPhase';
        end
        function Flash(obj, Status)
            % Virtual call            
        end       
        function GetPhaseEnthalpies(obj, Status)
            % due to including capillary pressure functions, P_2 is initialized
            p = Status.Properties('P_2').Value;
            for i=1:obj.NofPhases
                PhaseEnthalpy = Status.Properties(['h_',num2str(i)]);
                PhaseEnthalpy.Value = obj.Phases(i).ComputePhaseEnthalpy(i, p);    
            end
        end
        function GetPhaseIndex(obj, Status)
            h = Status.Properties('hTfluid').Value;
            PhaseIndex = Status.Properties('PhaseStatus');
            temp = zeros(length(h),1);
            temp(h <= Status.Properties('h_1').Value) = 1; % in liquid phase
            temp(h > Status.Properties('h_1').Value & h < Status.Properties('h_2').Value) = 2; % in two phase
            temp(h >= Status.Properties('h_2').Value) = 3; % in vapor phase
            PhaseIndex.Value = temp;
        end         
        function GetTemperature(obj, Status)
            PhaseIndex = Status.Properties('PhaseStatus').Value;
            p = Status.Properties('P_2').Value;
            h = Status.Properties('hTfluid').Value;
            % correct h for phase enthalpy of water in two-phase region
            h(PhaseIndex == 2) = Status.Properties('h_1').Value(PhaseIndex == 2);          
            T = Status.Properties('T');
            T.Value = obj.Phases(1).ComputeTemperature(PhaseIndex, p, h);
        end
        function GetPhaseViscosities(obj, Status)
            PhaseIndex = Status.Properties('PhaseStatus').Value;
            T = Status.Properties('T').Value;
            for i=1:obj.NofPhases
                mu = Status.Properties(['mu_',num2str(i)]);
                mu.Value = obj.Phases(i).ComputeViscosity(i, PhaseIndex, T);
            end
        end
        function GetPhaseDensities(obj, Status)
            PhaseIndex = Status.Properties('PhaseStatus').Value;
            p = Status.Properties('P_2').Value;
            h = Status.Properties('hTfluid').Value;            
            % compute saturation pressure
            psat = obj.Phases(1).ComputeSaturationPressure(Status, PhaseIndex);
            % correct p,h for saturation p,h in case of two-phase
            p(PhaseIndex == 2) = psat;            
            % compute density for both phases
            for i=1:obj.NofPhases
                % coorect h for saturation enthalpy in case of two-phase
                h(PhaseIndex == 2) = Status.Properties(['h_',num2str(i)]).Value(PhaseIndex == 2);
                % compute density 
                rho = Status.Properties(['rho_', num2str(i)]);
                rho.Value = obj.Phases(i).ComputeDensities(i, PhaseIndex, p, h);
            end
        end        
        function ComputePhaseSaturations(obj, Status)
            PhaseIndex = Status.Properties('PhaseStatus').Value;
            h = Status.Properties('hTfluid').Value;
            S1 = Status.Properties('S_1');
            S2 = Status.Properties('S_2');            
            % saturation liquid phase
            S1.Value(PhaseIndex == 1) = 1;
            S1.Value(PhaseIndex == 2) = ( Status.Properties('rho_2').Value(PhaseIndex == 2).* ...
                                                ( Status.Properties('h_2').Value(PhaseIndex == 2) - h(PhaseIndex == 2)) ) ./ ...
                                                ( h(PhaseIndex == 2).*( Status.Properties('rho_1').Value(PhaseIndex == 2) - ...
                                                  Status.Properties('rho_2').Value(PhaseIndex == 2) ) - ...
                                                ( Status.Properties('h_1').Value(PhaseIndex == 2).* Status.Properties('rho_1').Value(PhaseIndex == 2) - ...
                                                  Status.Properties('h_2').Value(PhaseIndex == 2).* Status.Properties('rho_2').Value(PhaseIndex == 2) ) );
            S1.Value(PhaseIndex == 3) = 0;
            
            % saturation vapor phase
            S2.Value(PhaseIndex == 1) = 0;
            S2.Value(PhaseIndex == 2) = 1 - S1.Value(PhaseIndex == 2);
            S2.Value(PhaseIndex == 3) = 1;   
        end   
        function correctEnthalpy(obj, Status)
            PhaseIndex = Status.Properties('PhaseStatus').Value;
            MixtureEnthalpy = Status.Properties('hTfluid');
            % correct phase enthalpies per phase (phase enthalpy is correct for two-phase region)
            for i=1:obj.NofPhases
                PhaseEnthalpy = Status.Properties(['h_',num2str(i)]);
                if i == 1
                    PhaseEnthalpy.Value(PhaseIndex == 1) = MixtureEnthalpy.Value(PhaseIndex == 1);
                    PhaseEnthalpy.Value(PhaseIndex == 3) = 0;
                elseif i == 2
                    PhaseEnthalpy.Value(PhaseIndex == 1) = 0;
                    PhaseEnthalpy.Value(PhaseIndex == 3) = MixtureEnthalpy.Value(PhaseIndex == 3);
                end
            end
        end
        
        % Compute variables
        function ComputeThermalConductivity(obj, Medium)
            % (1-phi)*C_r + phi*Sw*C_w + phi*Ss*C_s : vector
            CondEff = Medium.State.Properties('CondEff');
            CondEff.Value = (1 - Medium.Por) .* Medium.K_Cond_rock;
            for i=1:obj.NofPhases
                cond = [0.35; 0.6];
                S = Medium.State.Properties(['S_',num2str(i)]);
                CondEff.Value = CondEff.Value + Medium.Por .* cond(i) .* S.Value;
            end
        end
        function ComputeRockEnthalpy(obj, Medium)
            Hr = Medium.State.Properties('hRock');
            Hr.Value = Medium.Cpr .* Medium.State.Properties('T').Value;            
        end
        function Mob = ComputePhaseMobilities(obj, Status)
            S1 = Status.Properties('S_1').Value; 
            Mob = zeros(length(S1), obj.NofPhases);
            kr = obj.RelPermModel.ComputeRelPerm(obj.Phases, S1);
            for i=1:obj.NofPhases
                mu = Status.Properties(['mu_',num2str(i)]);
                Mob(:,i) = kr(:,i)./mu.Value;
            end
            Mob(isnan(Mob))=0;
        end
        
        % Derivatives for Jacobian MB and EB
        function drhodp = ComputeDrhoDp(obj, Status)
            PhaseIndex = Status.Properties('PhaseStatus').Value;
            p = Status.Properties('P_2').Value;
            h = Status.Properties('hTfluid').Value;
            % compute saturation pressure
            psat = obj.Phases(1).ComputeSaturationPressure(Status, PhaseIndex);
            % correct p,h for saturation p,h in case of two-phase
            p(PhaseIndex == 2) = psat; 
            % compute derivative
            drhodp = zeros(length(p),obj.NofPhases); % change this to length(p) and length(h)
            for i=1:obj.NofPhases
                % correct h for saturation enthalpy in case of two-phase
                h(PhaseIndex == 2) = Status.Properties(['h_',num2str(i)]).Value(PhaseIndex == 2);
                drhodp(:,i) = obj.Phases(i).ComputeDrhoDp(i, PhaseIndex, h, p);
            end
        end                     % MB
        function drhodh = ComputeDrhoDh(obj, Status)
            PhaseIndex = Status.Properties('PhaseStatus').Value;
            p = Status.Properties('P_2').Value;
            h = Status.Properties('hTfluid').Value;
            % compute saturation pressure
            psat = obj.Phases(1).ComputeSaturationPressure(Status, PhaseIndex);
            % correct p,h for saturation p,h in case of two-phase
            p(PhaseIndex == 2) = psat; 
            % compute derivative
            drhodh = zeros(length(h),obj.NofPhases);
            for i=1:obj.NofPhases
                % correct h for saturation enthalpy in case of two-phase
                h(PhaseIndex == 2) = Status.Properties(['h_',num2str(i)]).Value(PhaseIndex == 2);
                drhodh(:,i) = obj.Phases(i).ComputeDrhoDh(i, PhaseIndex, p, h);
            end
        end                     % MB, EB  
        function dTdp = ComputeDTDp(obj, Status)
            PhaseIndex = Status.Properties('PhaseStatus').Value;
            p = Status.Properties('P_2').Value;
            h = Status.Properties('hTfluid').Value;
            dTdp = obj.Phases(1).ComputeDTDp(PhaseIndex, p, h);
        end
        function dTdh = ComputeDTDh(obj, Status)
            PhaseIndex = Status.Properties('PhaseStatus').Value;
            p = Status.Properties('P_2').Value;
            h = Status.Properties('hTfluid').Value;
            % correct h for water saturation enthalpy in case of two-phase
            h(PhaseIndex == 2) = Status.Properties('h_1').Value(PhaseIndex == 2);          
            dTdh = obj.Phases(1).ComputeDTDh(PhaseIndex, p, h);
        end
        function dSdh = ComputeDSDh(obj, Status)
            h = Status.Properties('hTfluid').Value;
            h_w = Status.Properties('h_1').Value;
            h_s = Status.Properties('h_2').Value;
            rho_w = Status.Properties('rho_1').Value;
            rho_s = Status.Properties('rho_2').Value;
            % initialize
            dSdh = zeros(length(h),obj.NofPhases);
            % compute derivative
            dSdh(:,1) = ( rho_w .* rho_s .* ( h_w - h_s ) ) ./ ...
                        ( h .* ( rho_w - rho_s ) - ( h_w .* rho_w - h_s .* rho_s ) ).^2;
            dSdh(:,2) = -1 .* dSdh(:,1);
        end
        function dSdp = ComputeDSDp(obj, Status)
            %
        end

        function dMobdp = ComputeDMobDp(obj,Status)
            p = Status.Properties('P_2').Value;
            h = Status.Properties('hTfluid').Value;
            dMobdp = zeros(length(p),obj.NofPhases);
%             dSdp = zeros(length(p),obj.NofPhases);
            
            S1 = Status.Properties('S_1').Value; 
            kr = obj.RelPermModel.ComputeRelPerm(obj.Phases, S1);
            for i=1:obj.NofPhases
%                 dSdp(:,i) = obj.Phases(i).ComputeDSDp(obj.Pgrid, obj.Hgrid, STable, h, p);
                
                dMobdp(:,i) = 0;     
%                 dMobdp(:,i) = ( dSdp(:,i) .* mu.Value - dmudp(:,i) .* kr(:,i) ) ./ mu.Value.^2;     
            end
            dMobdp(isnan(dMobdp))=0;
        end
        function dMobdh = ComputeDMobDh(obj,Status)
            p = Status.Properties('P_2').Value;
            h = Status.Properties('hTfluid').Value;
            dMobdh = zeros(length(h),obj.NofPhases);
            obj.dSdh = ComputeDSDh(obj, Status);
            
            S1 = Status.Properties('S_1').Value; 
            kr = obj.RelPermModel.ComputeRelPerm(obj.Phases, S1);
            for i=1:obj.NofPhases
                mu = Status.Properties(['mu_',num2str(i)]).Value;
%                 dSdh(:,i) = obj.Phases(i).ComputeDSDh(obj.Pgrid, obj.Hgrid, STable, h, p);

                dMobdh(:,i) = 0;
%                 dMobdh(:,i) = ( obj.dSdh(:,i) .* mu ) ./ mu.^2;     
            end  
            dMobdh(isnan(dMobdh))=0;
        end
%         function dMobdh = ComputeDMobDh(obj,Status)
%             dMobdh = zeros(length(obj.Hindex),obj.NofPhases);
% 
%             % Part 1
%             rho_2Table = obj.TablePH.('rho_2');
%             drho_2dh = obj.Phases(2).ComputeDrhoDh(obj.Pgrid, obj.Hgrid, rho_2Table, obj.h, obj.p);
%             H_2 = Status.Properties('h_2').Value;
%             rho_2 = Status.Properties('rho_2').Value;
%             Part1 = drho_2dh .* H_2 - drho_2dh .* obj.h - rho_2 .* 1;
%             
%             % Part 2
%             rho_1Table = obj.TablePH.('rho_1');
%             drho_1dh = obj.Phases(2).ComputeDrhoDh(obj.Pgrid, obj.Hgrid, rho_1Table, obj.h, obj.p);
%             rho_1 = Status.Properties('rho_1').Value;
%             Part2 = 1 .* rho_1 + obj.h .* drho_1dh - 1 .* rho_2 - obj.h .* drho_2dh;
%             
%             % Part 3
%             H_1 = Status.Properties('h_1').Value;
%             Part3 = H_1 .* drho_1dh - H_2 .* drho_2dh;
%             
%             % chain rule (or quotient rule)
%             A = Part1 .* ( obj.h .* (rho_1 - rho_2) - (H_1 .* rho_1 - H_2 .* rho_2) );
%             B = rho_2 .* (H_2 - obj.h) .* ( Part2 - Part3 );
%             C = ( obj.h .* (rho_1 - rho_2) - (H_1 .* rho_1 - H_2 .* rho_2) ).^2;
%             
%             dMobdh(:,1) = (A - B) ./ C;
%             
%             dMobdh(:,2) = -1 .* dMobdh(:,1);
%             
%             dMobdh(isnan(dMobdh))=0;
%         end
        
        
        % These depend on how we treat the conductive flux
        function d2Td2p = ComputeD2TD2p(obj, Status)
            %
        end
        function d2Td2h = ComputeD2TD2h(obj, Status)
            %
        end

%         % 2nd derivatives for inflexion point correction/detection
%         function d2rhod2h = ComputeD2rhoD2h(obj, Status)
%             p = Status.Properties('P_2').Value;
%             h = Status.Properties('hTfluid').Value;
%             d2rhod2h = zeros(length(obj.Hindex),obj.NofPhases);
%             for i=1:obj.NofPhases
%                 rhoTable = obj.TablePH.(['rho_', num2str(i)]);
%                 d2rhod2h(:,i) = obj.Phases(i).ComputeD2rhoD2h(obj.Pgrid, obj.Hgrid, rhoTable, h, p);
%             end
%         end
%         function d2Mobdh2 = ComputeD2MobDh2(obj, Status)
%             p = Status.Properties('P_2').Value;
%             h = Status.Properties('hTfluid').Value;
%             dmudh = zeros(length(obj.Hindex),obj.NofPhases);             
%             d2mudh2 = zeros(length(obj.Hindex),obj.NofPhases); 
%             d2Mobdh2 = zeros(length(obj.Hindex),obj.NofPhases);
%             S1 = Status.Properties('S_1').Value; 
%             kr = obj.RelPermModel.ComputeRelPerm(obj.Phases, S1);
%             for i=1:obj.NofPhases
%                 mu = Status.Properties(['mu_',num2str(i)]);
%                 muTable = obj.TablePH.(['mu_', num2str(i)]);
%                 dmudh(:,i) = obj.Phases(i).ComputeDmuDh(obj.Pgrid, obj.Hgrid, muTable, h, p);
%                 d2mudh2(:,i) = obj.Phases(i).ComputeD2muDh2(obj.Pgrid, obj.Hgrid, muTable, h, p);
%                 
%                 d2Mobdh2(:,i) = -1 .* kr(:,i) .* ( d2mudh2(:,i) ./ mu.Value.^2 - 2 .* dmudh(:,i) ./ mu.Value.^3 );
%             end  
%             d2Mobdh2(isnan(d2Mobdh2))=0;
%         end
        


        % Other       
        function InitializeInjectors(obj, Inj)
            % Assuming saturation of injection phase is 1.0. We are assuming we are only injecting 1 phase, i.e. water
            for i=1:length(Inj)
                Inj(i).z = 1;
%                 Inj(i).S = 1;
                Inj(i).Mob = zeros(length(Inj(i).p),obj.NofPhases);
                
                if strcmp(Inj(i).BC_Formulation, 'Temperature')
                    Inj(i).h(:, 1) = obj.Phases(1).ComputeWaterEnthalpy(Inj(i).p, Inj(i).T);
                    Inj(i).rho(:, 1)= obj.Phases(1).ComputeWaterDensity(Inj(i).p, Inj(i).T);
                    mu = obj.Phases(1).ComputeWaterViscosity(Inj(i).T);
                elseif strcmp(Inj(i).BC_Formulation, 'Enthalpy')
%                     Inj(i).T = obj.Phases(1).ComputeTemperature(PhaseIndex, Inj(i).p, Inj(i).h(:,1));
                    Inj(i).T = 273.15 - 2.41231 + 2.56222e-8.*(Inj(i).h(:,1).*1e4) + ...
                               -9.31415e-17.*(Inj(i).p.*1e1).^2 - 2.2568e-19.*(Inj(i).h(:,1).*1e4).^2;
                           
                    Inj(i).rho(:, 1) = ( 1.00207 + 4.42607e-11.*(Inj(i).p.*1e1) - 5.47456e-12.*(Inj(i).h(:,1).*1e4) + ...
                        5.02875e-21.*(Inj(i).h(:,1).*1e4).*(Inj(i).p.*1e1) - 1.24791e-21.*(Inj(i).h(:,1).*1e4).^2 ).*1e3;
                    mu = (241.4 .* 10.^(247.8./((Inj(i).T - 273.15) + 133.15)) ) .* 1e-4 .* 1e-3;
                end
                    
                Inj(i).h(:, 2:obj.NofPhases) = 0;
                Inj(i).rho(:, 2:obj.NofPhases) = 0; 
                Inj(i).Mob(:, 1) = 1/mu; % injecting only water means rel.perm of water is 1.0 
                Inj(i).Mob(:, 2:obj.NofPhases) = 0;            
            end
        end
        
        function v = ComputeVelocity(obj, Reservoir, mu)
%             virtual call
        end  
        
    end
end
