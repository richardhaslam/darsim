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
        Ptable = (1:1:220)'.*1e5; % The range of pressure for all the tables
        Htable = (20:10:4800)'*1e3; % The range of enthalpy for all the tables
        Pindex % The index of pressure value for lookup in the tables
        Hindex % The index of enthalpy value for lookup in the tables
        Pstepsize = 1e5;
        Hstepsize = 1e4; % implement this
        
        Pgrid
        Hgrid
        p
        h
        
    end
    methods
        function obj = Geothermal_Multiphase_fluid_model(n_phases)
            obj@fluid_model(n_phases, n_phases);
            obj.name = 'Geothermal_MultiPhase';
        end
        function Flash(obj, Status)
            % Virtual call            
        end
        
        function GetTableIndex(obj, Status)
            % due to including capillary pressure functions, P_2 is initialized
            [~,obj.Pindex] = ismember( round(Status.Properties('P_2').Value,-log10(obj.Pstepsize)), round(obj.Ptable,-log10(obj.Pstepsize)) );
            [~,obj.Hindex] = ismember( round(Status.Properties('hTfluid').Value,-log10(obj.Hstepsize)), round(obj.Htable,-log10(obj.Hstepsize)) );
        end  % to get indices for pressure and enthalpy solution values
        function ComputeTableGrid(obj)
            [obj.Hgrid,obj.Pgrid] = meshgrid(obj.Htable,obj.Ptable);
        end
        function GetPHValues(obj, Status)
            obj.p = Status.Properties('P_2').Value;
            obj.h = Status.Properties('hTfluid').Value;
%             obj.p = abs(Status.Properties('P_2').Value);
%             obj.h = abs(Status.Properties('hTfluid').Value);
        end
        
        % "obj.Pgrid, obj.Hgrid" should be reversed; for consistency only !!
        
        
        % Linear interpolation
        function GetPhaseDensities(obj, Status)
            for i=1:obj.NofPhases
                rho = Status.Properties(['rho_', num2str(i)]);
                rhoTable = obj.TablePH.(['rho_', num2str(i)]);
                rho.Value = obj.Phases(i).GetDensity(obj.Pgrid, obj.Hgrid, rhoTable, obj.h, obj.p);
            end
        end
        function GetPhaseSaturations(obj, Status)
            for i=1:obj.NofPhases
                S = Status.Properties(['S_',num2str(i)]); % S is just a pointer to a memory location here, so it has no value
                STable = obj.TablePH.(['S_',num2str(i)]); % |-> This is why you need num2str() again here
                S.Value = obj.Phases(i).GetSaturation(obj.Pgrid, obj.Hgrid, STable, obj.h, obj.p);
            end
        end
        function GetPhaseViscosities(obj, Status)
            for i=1:obj.NofPhases
                mu = Status.Properties(['mu_',num2str(i)]);
                muTable = obj.TablePH.(['mu_',num2str(i)]);
                mu.Value = obj.Phases(i).GetViscosity(obj.Pgrid, obj.Hgrid, muTable, obj.h, obj.p);
            end
        end
        function GetPhaseThermalConductivities(obj, Status) 
            for i=1:obj.NofPhases
                ThermCond = Status.Properties(['cond_',num2str(i)]);
                ThermCondTable = obj.TablePH.(['cond_',num2str(i)]);
                ThermCond.Value = obj.Phases(i).GetConductivity(obj.Pgrid, obj.Hgrid, ThermCondTable, obj.h, obj.p);
            end
        end
        function GetPhaseEnthalpies(obj, Status)
            MixtureEnthalpy = Status.Properties('hTfluid');
            for i=1:obj.NofPhases
                PhaseEnthalpy = Status.Properties(['h_',num2str(i)]);
                PhaseEnthalpyTable = obj.TablePH.(['H_',num2str(i)]);
                PhaseEnthalpy.Value = obj.Phases(i).GetPhaseEnthalpy(obj.Ptable, PhaseEnthalpyTable, obj.p);
                
                % Correcting the phase enthalpies in the single phase regions
                Index1 = find( Status.Properties(['S_',num2str(i)]).Value == 1 );
                Index2 = find( Status.Properties(['S_',num2str(i)]).Value == 0 );
                PhaseEnthalpy.Value(Index1) = MixtureEnthalpy.Value(Index1);
                PhaseEnthalpy.Value(Index2) = 0;
            end
        end

        % Get Total/General properties
        function GetTemperature(obj, Status)
            T = Status.Properties('T');
            T.Value = interp2(obj.Hgrid, obj.Pgrid, obj.TablePH.Temperature, obj.h, obj.p, 'linear');
        end
        function GetTotalDensity(obj, Status)
            rhoT = Status.Properties('rhoT');
            rhoT.Value = interp2(obj.Hgrid, obj.Pgrid, obj.TablePH.rhoT, obj.h, obj.p, 'linear');
        end
        
        % Compute variables
        function ComputeThermalConductivity(obj, Reservoir)
            % (1-phi)*C_r + phi*Sw*C_w + phi*Ss*C_s : vector
            CondEff = Reservoir.State.Properties('CondEff');
            CondEff.Value = (1 - Reservoir.Por) .* Reservoir.K_Cond_rock;
            for i=1:obj.NofPhases
                cond = Reservoir.State.Properties(['cond_',num2str(i)]);
                S = Reservoir.State.Properties(['S_',num2str(i)]);
                CondEff.Value = CondEff.Value + Reservoir.Por .* cond.Value .* S.Value;
            end
        end
        function ComputeRockEnthalpy(obj, Reservoir)
            Hr = Reservoir.State.Properties('hRock');
            Hr.Value = Reservoir.Cpr .* Reservoir.State.Properties('T').Value;            
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
        function drhodp = ComputeDrhoDp(obj)
            drhodp = zeros(length(obj.Pindex),obj.NofPhases); % change this to length(p) and length(h)
            for i=1:obj.NofPhases
                rhoTable = obj.TablePH.(['rho_', num2str(i)]);
                drhodp(:,i) = obj.Phases(i).ComputeDrhoDp(obj.Pgrid, obj.Hgrid, rhoTable, obj.h, obj.p);
            end
        end                     % MB
        function drhodh = ComputeDrhoDh(obj)
            drhodh = zeros(length(obj.Hindex),obj.NofPhases);
            for i=1:obj.NofPhases
                rhoTable = obj.TablePH.(['rho_', num2str(i)]);
                drhodh(:,i) = obj.Phases(i).ComputeDrhoDh(obj.Pgrid, obj.Hgrid, rhoTable, obj.h, obj.p);
            end
        end                     % MB, EB
        function drhoTdp = ComputeDrhoTDp(obj) 
            [~,table_drhoTdp] = gradient(obj.TablePH.rhoT,obj.Pstepsize); 
            drhoTdp = interp2(obj.Hgrid, obj.Pgrid, table_drhoTdp, obj.h, obj.p, 'linear');     
        end                  % MB
        function drhoTdh = ComputeDrhoTDh(obj) 
            [table_drhoTdh,~] = gradient(obj.TablePH.rhoT,obj.Hstepsize); 
            drhoTdh = interp2(obj.Hgrid, obj.Pgrid, table_drhoTdh, obj.h, obj.p, 'linear');  
        end                  % MB
        function drho_times_Sdp = ComputeDrho_times_SDp(obj) 
            drho_times_Sdp = zeros(length(obj.Pindex),obj.NofPhases);
            for i=1:obj.NofPhases
                rho_times_STable = obj.TablePH.(['rho_times_S_',num2str(i)]);
                drho_times_Sdp(:,i) = obj.Phases(i).ComputeDrho_times_SDp(obj.Pgrid, obj.Hgrid, rho_times_STable, obj.h, obj.p);
            end
        end    % MB, EB
        function drho_times_Sdh = ComputeDrho_times_SDh(obj)
            drho_times_Sdh = zeros(length(obj.Hindex),obj.NofPhases);
            for i=1:obj.NofPhases
                rho_times_STable = obj.TablePH.(['rho_times_S_',num2str(i)]);
                drho_times_Sdh(:,i) = obj.Phases(i).ComputeDrho_times_SDh(obj.Pgrid, obj.Hgrid, rho_times_STable, obj.h, obj.p);
            end            
        end     % MB, EB
        function drho_times_hdp = ComputeDrho_times_hDp(obj)
            drho_times_hdp = zeros(length(obj.Pindex),obj.NofPhases);
            for i=1:obj.NofPhases
                rho_times_hTable = obj.TablePH.(['rho_times_H_',num2str(i)]);
                drho_times_hdp(:,i) = obj.Phases(i).ComputeDrho_times_hDp(obj.Pgrid, obj.Hgrid, rho_times_hTable, obj.h, obj.p);
            end
        end     % EB     
        function drho_times_hdh = ComputeDrho_times_hDh(obj)
            drho_times_hdh = zeros(length(obj.Hindex),obj.NofPhases);
            for i=1:obj.NofPhases
                rho_times_hTable = obj.TablePH.(['rho_times_H_',num2str(i)]);
                drho_times_hdh(:,i) = obj.Phases(i).ComputeDrho_times_hDh(obj.Pgrid, obj.Hgrid, rho_times_hTable, obj.h, obj.p);
            end
        end     % EB     
        function dMobdp = ComputeDMobDp(obj,Status)
            dmudp = zeros(length(obj.Pindex),obj.NofPhases); 
            dMobdp = zeros(length(obj.Pindex),obj.NofPhases);
            dSdp = zeros(length(obj.Pindex),obj.NofPhases);
            
            S1 = Status.Properties('S_1').Value; 
            kr = obj.RelPermModel.ComputeRelPerm(obj.Phases, S1);
            for i=1:obj.NofPhases
                mu = Status.Properties(['mu_',num2str(i)]);
                muTable = obj.TablePH.(['mu_', num2str(i)]);
                STable = obj.TablePH.(['S_', num2str(i)]);
                
                dmudp(:,i) = obj.Phases(i).ComputeDmuDp(obj.Pgrid, obj.Hgrid, muTable, obj.h, obj.p);
                dSdp(:,i) = obj.Phases(i).ComputeDSDp(obj.Pgrid, obj.Hgrid, STable, obj.h, obj.p);
                
                dMobdp(:,i) = ( -1 .* dmudp(:,i) .* kr(:,i) ) ./ mu.Value.^2;     
%                 dMobdp(:,i) = ( dSdp(:,i) .* mu.Value - dmudp(:,i) .* kr(:,i) ) ./ mu.Value.^2;     
            end
            dMobdp(isnan(dMobdp))=0;
        end
        function dMobdh = ComputeDMobDh(obj,Status)
            dmudh = zeros(length(obj.Hindex),obj.NofPhases); 
            dMobdh = zeros(length(obj.Hindex),obj.NofPhases);
            dSdh = zeros(length(obj.Pindex),obj.NofPhases);

            S1 = Status.Properties('S_1').Value; 
            kr = obj.RelPermModel.ComputeRelPerm(obj.Phases, S1);
            for i=1:obj.NofPhases
                mu = Status.Properties(['mu_',num2str(i)]);
                muTable = obj.TablePH.(['mu_', num2str(i)]);
                STable = obj.TablePH.(['S_', num2str(i)]);

                dmudh(:,i) = obj.Phases(i).ComputeDmuDh(obj.Pgrid, obj.Hgrid, muTable, obj.h, obj.p);
                dSdh(:,i) = obj.Phases(i).ComputeDSDh(obj.Pgrid, obj.Hgrid, STable, obj.h, obj.p);

%                 dMobdh(:,i) = ( -1 .* dmudh(:,i) .* kr(:,i) ) ./ mu.Value.^2;
                dMobdh(:,i) = ( dSdh(:,i) .* mu.Value - dmudh(:,i) .* kr(:,i) ) ./ mu.Value.^2;     
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
        function dSdp = ComputeDSDp(obj)
            dSdp = zeros(length(obj.Pindex),obj.NofPhases);
            for i=1:obj.NofPhases
                STable = obj.TablePH.(['S_', num2str(i)]);
                dSdp(:,i) = obj.Phases(i).ComputeDSDp(obj.Pgrid, obj.Hgrid, STable, obj.h, obj.p);
            end
        end
        function dSdh = ComputeDSDh(obj)
            dSdh = zeros(length(obj.Hindex),obj.NofPhases);
            for i=1:obj.NofPhases
                STable = obj.TablePH.(['S_', num2str(i)]);
                dSdh(:,i) = obj.Phases(i).ComputeDSDh(obj.Pgrid, obj.Hgrid, STable, obj.h, obj.p);
            end
        end
        function drhoHSdp = ComputeDrhoHSDp(obj)
            drhoHSdp = zeros(length(obj.Pindex),obj.NofPhases);
            for i=1:obj.NofPhases
                rhoHSTable = obj.TablePH.(['rhoHS_',num2str(i)]);
                drhoHSdp(:,i) = obj.Phases(i).ComputeDrhoHSDp(obj.Pgrid, obj.Hgrid, rhoHSTable, obj.h, obj.p);
            end
        end
        function drhoHSdh = ComputeDrhoHSDh(obj)
            drhoHSdh = zeros(length(obj.Hindex),obj.NofPhases);
            for i=1:obj.NofPhases
                rhoHSTable = obj.TablePH.(['rhoHS_',num2str(i)]);
                drhoHSdh(:,i) = obj.Phases(i).ComputeDrhoHSDh(obj.Pgrid, obj.Hgrid, rhoHSTable, obj.h, obj.p);
            end            
        end
        function dTdh = ComputeDTDh(obj)
            TTable = obj.TablePH.Temperature;
            dTdh = obj.Phases(1).ComputeDTDh(obj.Pgrid, obj.Hgrid, TTable, obj.h, obj.p);
        end
        function dTdp = ComputeDTDp(obj)
            TTable = obj.TablePH.Temperature; % What if we just pass obj.TablePH.Temperature to the function ??
            dTdp = obj.Phases(1).ComputeDTDp(obj.Pgrid, obj.Hgrid, TTable, obj.h, obj.p);
        end
        
        % These depend on how we treat the conductive flux
        function d2Td2p = ComputeD2TD2p(obj)
            TTable = obj.TablePH.Temperature; % What if we just pass obj.TablePH.Temperature to the function ??
            d2Td2p = obj.Phases(1).ComputeD2TD2p(obj.Pgrid, obj.Hgrid, TTable, obj.h, obj.p);
        end
        function d2Td2h = ComputeD2TD2h(obj)
            TTable = obj.TablePH.Temperature;
            d2Td2h = obj.Phases(1).ComputeD2TD2h(obj.Pgrid, obj.Hgrid, TTable, obj.h, obj.p);
        end

        % 2nd derivatives for inflexion point correction/detection
        function d2rhod2h = ComputeD2rhoD2h(obj)
            d2rhod2h = zeros(length(obj.Hindex),obj.NofPhases);
            for i=1:obj.NofPhases
                rhoTable = obj.TablePH.(['rho_', num2str(i)]);
                d2rhod2h(:,i) = obj.Phases(i).ComputeD2rhoD2h(obj.Pgrid, obj.Hgrid, rhoTable, obj.h, obj.p);
            end
        end
        function d2Mobdh2 = ComputeD2MobDh2(obj)
            dmudh = zeros(length(obj.Hindex),obj.NofPhases);             
            d2mudh2 = zeros(length(obj.Hindex),obj.NofPhases); 
            d2Mobdh2 = zeros(length(obj.Hindex),obj.NofPhases);
            S1 = Status.Properties('S_1').Value; 
            kr = obj.RelPermModel.ComputeRelPerm(obj.Phases, S1);
            for i=1:obj.NofPhases
                mu = Status.Properties(['mu_',num2str(i)]);
                muTable = obj.TablePH.(['mu_', num2str(i)]);
                dmudh(:,i) = obj.Phases(i).ComputeDmuDh(obj.Pgrid, obj.Hgrid, muTable, obj.h, obj.p);
                d2mudh2(:,i) = obj.Phases(i).ComputeD2muDh2(obj.Pgrid, obj.Hgrid, muTable, obj.h, obj.p);
                
                d2Mobdh2(:,i) = -1 .* kr(:,i) .* ( d2mudh2(:,i) ./ mu.Value.^2 - 2 .* dmudh(:,i) ./ mu.Value.^3 );
            end  
            d2Mobdh2(isnan(d2Mobdh2))=0;
        end
        
        % Other       
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
        
        function v = ComputeVelocity(obj, Reservoir, mu)
%             virtual call
        end  
        
    end
end
