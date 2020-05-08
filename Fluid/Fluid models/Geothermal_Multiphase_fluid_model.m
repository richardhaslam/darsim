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
%         Ptable = (1:0.1:220)'.*1e5; % The range of pressure for all the tables
%         Htable = (20:1:4800)'*1e3; % The range of enthalpy for all the tables
        Ptable = (1:0.1:220)'.*1e5; % The range of pressure for all the tables
        Htable = (20:1:4800)'*1e3; % The range of enthalpy for all the tables
        Pindex % The index of pressure value for lookup in the tables
        Hindex % The index of enthalpy value for lookup in the tables
%         Pstepsize = 1e4;
%         Hstepsize = 1e3; % implement this
        Pstepsize = 1e4;
        Hstepsize = 1e3; % implement this
        
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

%         function GetTableIndex(obj, Status)
%             % due to including capillary pressure functions, P_2 is initialized
%             [~,obj.Pindex] = ismember( round(Status.Properties('P_2').Value,-log10(obj.Pstepsize)), round(obj.Ptable,-log10(obj.Pstepsize)) );
%             [~,obj.Hindex] = ismember( round(Status.Properties('hTfluid').Value,-log10(obj.Hstepsize)), round(obj.Htable,-log10(obj.Hstepsize)) );
%         end  % to get indices for pressure and enthalpy solution values

        % Get Phase properties
% %         function GetPhaseDensities(obj, Status) 
% %             for i=1:obj.NofPhases
% %                 rho = Status.Properties(['rho_', num2str(i)]);
% %                 rhoTable = obj.TablePH.(['rho_', num2str(i)]);
% %                 rho.Value = obj.Phases(i).GetDensity(obj.Pindex, obj.Hindex, rhoTable);
% %             end
% %         end        
% %         function GetPhaseSaturations(obj, Status)
% %             for i=1:obj.NofPhases
% %                 S = Status.Properties(['S_',num2str(i)]); % S is just a pointer to a memory location here, so it has no value
% %                 STable = obj.TablePH.(['S_',num2str(i)]); % |-> This is why you need num2str() again here
% %                 S.Value = obj.Phases(i).GetSaturation(obj.Pindex, obj.Hindex, STable); 
% %             end
% %         end
% %         function GetPhaseViscosities(obj, Status)
% %             for i=1:obj.NofPhases
% %                 mu = Status.Properties(['mu_',num2str(i)]);
% %                 muTable = obj.TablePH.(['mu_',num2str(i)]);
% %                 mu.Value = obj.Phases(i).GetViscosity(obj.Pindex, obj.Hindex, muTable);
% %             end
% %         end
% %         function GetPhaseThermalConductivities(obj, Status) 
% %             for i=1:obj.NofPhases
% %                 ThermCond = Status.Properties(['cond_',num2str(i)]);
% %                 ThermCondTable = obj.TablePH.(['cond_',num2str(i)]);
% %                 ThermCond.Value = obj.Phases(i).GetConductivity(obj.Pindex, obj.Hindex, ThermCondTable);
% %             end
% %         end
% %         function GetPhaseEnthalpies(obj, Status)
% %             for i=1:obj.NofPhases
% %                 PhaseEnthalpy = Status.Properties(['h_',num2str(i)]);
% %                 PhaseEnthalpyTable = obj.TablePH.(['H_',num2str(i)]);
% %                 PhaseEnthalpy.Value = obj.Phases(i).GetPhaseEnthalpy(obj.Pindex, PhaseEnthalpyTable);    
% %             end
% %         end
% %     
% %         % Get Total/General properties
% %         function GetTemperature(obj, Status)
% %             T = Status.Properties('T');
% %             T.Value = obj.TablePH.Temperature(sub2ind(size(obj.TablePH.Temperature), obj.Pindex, obj.Hindex));
% %         end
% %         function GetTotalDensity(obj, Status)
% %             rhoT = Status.Properties('rhoT');
% %             rhoT.Value = obj.TablePH.rhoT( sub2ind(size(obj.TablePH.rhoT), obj.Pindex, obj.Hindex) ); 
% %         end
% %         
% %         % Compute variables
% %         function ComputeThermalConductivity(obj, Reservoir)
% %             % (1-phi)*C_r + phi*Sw*C_w + phi*Ss*C_s : vector
% %             CondEff = Reservoir.State.Properties('CondEff');
% %             CondEff.Value = (1 - Reservoir.Por) .* Reservoir.K_Cond_rock;
% %             for i=1:obj.NofPhases
% %                 cond = Reservoir.State.Properties(['cond_',num2str(i)]);
% %                 S = Reservoir.State.Properties(['S_',num2str(i)]);
% %                 CondEff.Value = CondEff.Value + Reservoir.Por .* cond.Value .* S.Value;
% %             end
% %         end
% %         function ComputeRockEnthalpy(obj, Reservoir)
% %             Hr = Reservoir.State.Properties('hRock');
% %             Hr.Value = Reservoir.Cpr .* Reservoir.State.Properties('T').Value;            
% %         end
% %         function Mob = ComputePhaseMobilities(obj, Status)
% %             S1 = Status.Properties('S_1').Value; 
% %             Mob = zeros(length(S1), obj.NofPhases);
% %             kr = obj.RelPermModel.ComputeRelPerm(obj.Phases, S1);
% %             for i=1:obj.NofPhases
% %                 mu = Status.Properties(['mu_',num2str(i)]);
% %                 Mob(:,i) = kr(:,i)./mu.Value;
% %             end
% %             Mob(isnan(Mob))=0;
% %         end
% %         
% %         % Derivatives for Jacobian MB and EB
% %         function drhodp = ComputeDrhoDp(obj)
% %             drhodp = zeros(length(obj.Pindex),obj.NofPhases);
% %             for i=1:obj.NofPhases
% %                 rhoTable = obj.TablePH.(['rho_', num2str(i)]);
% %                 drhodp(:,i) = obj.Phases(i).ComputeDrhoDp(obj.Pindex, obj.Hindex, rhoTable);
% %             end
% %         end                     % MB
% %         function drhodh = ComputeDrhoDh(obj)
% %             drhodh = zeros(length(obj.Hindex),obj.NofPhases);
% %             for i=1:obj.NofPhases
% %                 rhoTable = obj.TablePH.(['rho_', num2str(i)]);
% %                 drhodh(:,i) = obj.Phases(i).ComputeDrhoDh(obj.Pindex, obj.Hindex, rhoTable);
% %             end
% %         end                     % MB, EB
% %         function drhoTdp = ComputeDrhoTDp(obj) 
% %             [~,table_drhoTdp] = gradient(obj.TablePH.rhoT,obj.Pstepsize); 
% %             drhoTdp = table_drhoTdp(sub2ind(size(table_drhoTdp), obj.Pindex, obj.Hindex));      
% %         end                  % MB
% %         function drhoTdh = ComputeDrhoTDh(obj) 
% %             [table_drhoTdh,~] = gradient(obj.TablePH.rhoT,obj.Hstepsize); 
% %             drhoTdh = table_drhoTdh(sub2ind(size(table_drhoTdh), obj.Pindex, obj.Hindex));
% %         end                  % MB
% %         function drho_times_Sdp = ComputeDrho_times_SDp(obj) 
% %             drho_times_Sdp = zeros(length(obj.Pindex),obj.NofPhases);
% %             for i=1:obj.NofPhases
% %                 rho_times_STable = obj.TablePH.(['rho_times_S_',num2str(i)]);
% %                 drho_times_Sdp(:,i) = obj.Phases(i).ComputeDrho_times_SDp(obj.Pindex, obj.Hindex, rho_times_STable);
% %             end
% %         end    % MB, EB
% %         function drho_times_Sdh = ComputeDrho_times_SDh(obj)
% %             drho_times_Sdh = zeros(length(obj.Hindex),obj.NofPhases);
% %             for i=1:obj.NofPhases
% %                 rho_times_STable = obj.TablePH.(['rho_times_S_',num2str(i)]);
% %                 drho_times_Sdh(:,i) = obj.Phases(i).ComputeDrho_times_SDh(obj.Pindex, obj.Hindex, rho_times_STable);
% %             end            
% %         end     % MB, EB
% %         function drho_times_hdp = ComputeDrho_times_hDp(obj)
% %             drho_times_hdp = zeros(length(obj.Pindex),obj.NofPhases);
% %             for i=1:obj.NofPhases
% %                 rho_times_hTable = obj.TablePH.(['rho_times_H_',num2str(i)]);
% %                 drho_times_hdp(:,i) = obj.Phases(i).ComputeDrho_times_hDp(obj.Pindex, obj.Hindex, rho_times_hTable);
% %             end
% %         end     % EB     
% %         function drho_times_hdh = ComputeDrho_times_hDh(obj)
% %             drho_times_hdh = zeros(length(obj.Hindex),obj.NofPhases);
% %             for i=1:obj.NofPhases
% %                 rho_times_hTable = obj.TablePH.(['rho_times_H_',num2str(i)]);
% %                 drho_times_hdh(:,i) = obj.Phases(i).ComputeDrho_times_hDh(obj.Pindex, obj.Hindex, rho_times_hTable);
% %             end
% %         end     % EB     
% %         function dMobdp = ComputeDMobDp(obj,Status)
% %             dmudp = zeros(length(obj.Pindex),obj.NofPhases); 
% %             dMobdp = zeros(length(obj.Pindex),obj.NofPhases);
% %             S1 = Status.Properties('S_1').Value; 
% %             kr = obj.RelPermModel.ComputeRelPerm(obj.Phases, S1);
% %             for i=1:obj.NofPhases
% %                 mu = Status.Properties(['mu_',num2str(i)]);
% %                 muTable = obj.TablePH.(['mu_', num2str(i)]);
% %                 dmudp(:,i) = obj.Phases(i).ComputeDmuDp(obj.Pindex, obj.Hindex, muTable);
% %                 dMobdp(:,i) = ( -1 .* dmudp(:,i) .* kr(:,i) ) ./ mu.Value.^2;     
% %             end
% %             dMobdp(isnan(dMobdp))=0;
% %         end
% %         function dMobdh = ComputeDMobDh(obj,Status)
% %             dmudh = zeros(length(obj.Hindex),obj.NofPhases); 
% %             dMobdh = zeros(length(obj.Hindex),obj.NofPhases);
% %             S1 = Status.Properties('S_1').Value; 
% %             kr = obj.RelPermModel.ComputeRelPerm(obj.Phases, S1);
% %             for i=1:obj.NofPhases
% %                 mu = Status.Properties(['mu_',num2str(i)]);
% %                 muTable = obj.TablePH.(['mu_', num2str(i)]);
% %                 dmudh(:,i) = obj.Phases(i).ComputeDmuDh(obj.Pindex, obj.Hindex, muTable);
% %                 dMobdh(:,i) = ( -1 .* dmudh(:,i) .* kr(:,i) ) ./ mu.Value.^2;
% %             end  
% %             dMobdh(isnan(dMobdh))=0;
% %         end
% %         function dSdp = ComputeDSDp(obj)
% %             dSdp = zeros(length(obj.Pindex),obj.NofPhases);
% %             for i=1:obj.NofPhases
% %                 STable = obj.TablePH.(['S_', num2str(i)]);
% %                 dSdp(:,i) = obj.Phases(i).ComputeDSDp(obj.Pindex, obj.Hindex, STable);
% %             end
% %         end
% %         function dSdh = ComputeDSDh(obj)
% %             dSdh = zeros(length(obj.Hindex),obj.NofPhases);
% %             for i=1:obj.NofPhases
% %                 STable = obj.TablePH.(['S_', num2str(i)]);
% %                 dSdh(:,i) = obj.Phases(i).ComputeDSDh(obj.Pindex, obj.Hindex, STable);
% %             end
% %         end
% %         function drhoHSdp = ComputeDrhoHSDp(obj)
% %             drhoHSdp = zeros(length(obj.Pindex),obj.NofPhases);
% %             for i=1:obj.NofPhases
% %                 rhoHSTable = obj.TablePH.(['rhoHS_',num2str(i)]);
% %                 drhoHSdp(:,i) = obj.Phases(i).ComputeDrhoHSDp(obj.Pindex, obj.Hindex, rhoHSTable);
% %             end
% %         end
% %         function drhoHSdh = ComputeDrhoHSDh(obj)
% %             drhoHSdh = zeros(length(obj.Hindex),obj.NofPhases);
% %             for i=1:obj.NofPhases
% %                 rhoHSTable = obj.TablePH.(['rhoHS_',num2str(i)]);
% %                 drhoHSdh(:,i) = obj.Phases(i).ComputeDrhoHSDh(obj.Pindex, obj.Hindex, rhoHSTable);
% %             end            
% %         end
% %         function dTdh = ComputeDTDh(obj)
% %             TTable = obj.TablePH.Temperature;
% %             dTdh = obj.Phases(1).ComputeDTDh(obj.Pindex, obj.Hindex, TTable);
% %         end
% %         function dTdp = ComputeDTDp(obj)
% %             TTable = obj.TablePH.Temperature; % What if we just pass obj.TablePH.Temperature to the function ??
% %             dTdp = obj.Phases(1).ComputeDTDp(obj.Pindex, obj.Hindex, TTable);
% %         end
% % 
% %         % These depend on how we treat the conductive flux
% %         function d2Td2p = ComputeD2TD2p(obj)
% %             TTable = obj.TablePH.Temperature; % What if we just pass obj.TablePH.Temperature to the function ??
% %             d2Td2p = obj.Phases(1).ComputeD2TD2p(obj.Pindex, obj.Hindex, TTable);
% %         end
% %         function d2Td2h = ComputeD2TD2h(obj)
% %             TTable = obj.TablePH.Temperature;
% %             d2Td2h = obj.Phases(1).ComputeD2TD2h(obj.Pindex, obj.Hindex, TTable);
% %         end
% % 
% %         % 2nd derivatives for inflexion point correction/detection
% %         function d2rhod2h = ComputeD2rhoD2h(obj)
% %             d2rhod2h = zeros(length(obj.Hindex),obj.NofPhases);
% %             for i=1:obj.NofPhases
% %                 rhoTable = obj.TablePH.(['rho_', num2str(i)]);
% %                 d2rhod2h(:,i) = obj.Phases(i).ComputeD2rhoD2h(obj.Pindex, obj.Hindex, rhoTable);
% %             end
% %         end
% %         function d2Mobdh2 = ComputeD2MobDh2(obj)
% %             dmudh = zeros(length(obj.Hindex),obj.NofPhases);             
% %             d2mudh2 = zeros(length(obj.Hindex),obj.NofPhases); 
% %             d2Mobdh2 = zeros(length(obj.Hindex),obj.NofPhases);
% %             S1 = Status.Properties('S_1').Value; 
% %             kr = obj.RelPermModel.ComputeRelPerm(obj.Phases, S1);
% %             for i=1:obj.NofPhases
% %                 mu = Status.Properties(['mu_',num2str(i)]);
% %                 muTable = obj.TablePH.(['mu_', num2str(i)]);
% %                 dmudh(:,i) = obj.Phases(i).ComputeDmuDh(obj.Pindex, obj.Hindex, muTable);
% %                 d2mudh2(:,i) = obj.Phases(i).ComputeD2muDh2(obj.Pindex, obj.Hindex, muTable);
% %                 
% %                 d2Mobdh2(:,i) = -1 .* kr(:,i) .* ( d2mudh2(:,i) ./ mu.Value.^2 - 2 .* dmudh(:,i) ./ mu.Value.^3 );
% %             end  
% %             d2Mobdh2(isnan(d2Mobdh2))=0;
% %         end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% LINEAR INTERPOLATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
            T.Value = interp2(obj.Hgrid, obj.Pgrid, obj.TablePH.Temperature, obj.h, obj.p);
        end
        function GetTotalDensity(obj, Status)
            rhoT = Status.Properties('rhoT');
            rhoT.Value = interp2(obj.Hgrid, obj.Pgrid, obj.TablePH.rhoT, obj.h, obj.p);
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
            drhodp = zeros(length(obj.Pindex),obj.NofPhases);
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
            drhoTdp = interp2(obj.Hgrid, obj.Pgrid, table_drhoTdp, obj.h, obj.p);     
        end                  % MB
        function drhoTdh = ComputeDrhoTDh(obj) 
            [table_drhoTdh,~] = gradient(obj.TablePH.rhoT,obj.Hstepsize); 
            drhoTdh = interp2(obj.Hgrid, obj.Pgrid, table_drhoTdh, obj.h, obj.p);  
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
            S1 = Status.Properties('S_1').Value; 
            kr = obj.RelPermModel.ComputeRelPerm(obj.Phases, S1);
            for i=1:obj.NofPhases
                mu = Status.Properties(['mu_',num2str(i)]);
                muTable = obj.TablePH.(['mu_', num2str(i)]);
                dmudp(:,i) = obj.Phases(i).ComputeDmuDp(obj.Pgrid, obj.Hgrid, muTable, obj.h, obj.p);
                dMobdp(:,i) = ( -1 .* dmudp(:,i) .* kr(:,i) ) ./ mu.Value.^2;     
            end
            dMobdp(isnan(dMobdp))=0;
        end
        function dMobdh = ComputeDMobDh(obj,Status)
            dmudh = zeros(length(obj.Hindex),obj.NofPhases); 
            dMobdh = zeros(length(obj.Hindex),obj.NofPhases);
            S1 = Status.Properties('S_1').Value; 
            kr = obj.RelPermModel.ComputeRelPerm(obj.Phases, S1);
            for i=1:obj.NofPhases
                mu = Status.Properties(['mu_',num2str(i)]);
                muTable = obj.TablePH.(['mu_', num2str(i)]);
                dmudh(:,i) = obj.Phases(i).ComputeDmuDh(obj.Pgrid, obj.Hgrid, muTable, obj.h, obj.p);
                dMobdh(:,i) = ( -1 .* dmudh(:,i) .* kr(:,i) ) ./ mu.Value.^2;
            end  
            dMobdh(isnan(dMobdh))=0;
        end
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
                Inj(i).x = [1 0];
                Inj(i).S = 1;
                Inj(i).Mob = zeros(length(Inj(i).p),obj.NofPhases);
                
                if strcmp(Inj(i).BC_Formulation, 'Temperature')
                    Inj(i).h(:, 1) = obj.Phases(1).ComputeWaterEnthalpy(Inj(i).p, Inj(i).T);
                    Inj(i).rho(:, 1)= obj.Phases(1).ComputeWaterDensity(Inj(i).p, Inj(i).T);
                    mu = obj.Phases(1).ComputeWaterViscosity(Inj(i).T);
                elseif strcmp(Inj(i).BC_Formulation, 'Enthalpy')
                    Inj(i).T = interp2(obj.Hgrid, obj.Pgrid, obj.TablePH.Temperature, Inj(i).h(:,1), Inj(i).p);
                    Inj(i).rho(:, 1) = obj.Phases(1).GetDensity(obj.Pgrid, obj.Hgrid, obj.TablePH.rho_1, Inj(i).h(:,1), Inj(i).p);
                    mu = obj.Phases(1).GetViscosity(obj.Pgrid, obj.Hgrid, obj.TablePH.mu_1, Inj(i).h(:,1), Inj(i).p);
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
        
        
        
        
        %%% ORIGINALLY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function GetPhaseInternalEnergies(obj, Status)
            for i=1:obj.NofPhases
                U = Status.Properties(['U_',num2str(i)]);
                UTable = obj.TablePH.(['U_',num2str(i)]);
                U.Value = obj.Phases(i).GetInternalEnergy(obj.Pindex, obj.Hindex, UTable);
            end
        end

        % Derivatives for Thermal Conductivity tensor
        function ds_times_conddp = ComputeDs_times_condDp(obj)
            ds_times_conddp = zeros(length(obj.Pindex),obj.NofPhases);
            for i=1:NofPhases
                s_times_condTable = obj.TablePH.(['s_times_cond_',num2str(i)]);
                ds_times_conddp(:,i) = obj.Phases(i).ComputeDs_times_condDp(obj.Pindex, obj.Hindex, s_times_condTable);
            end            
        end
        function ds_times_conddh = ComputeDs_times_condDh(obj)
            ds_times_conddh = zeros(length(obj.Hindex),obj.NofPhases);
            for i=1:NofPhases
                s_times_condTable = obj.TablePH.(['s_times_cond_',num2str(i)]);
                ds_times_conddh(:,i) = obj.Phases(i).ComputeDs_times_condDh(obj.Pindex, obj.Hindex, s_times_condTable);
            end
        end
        function dDdp = ComputeDDDp(obj, Status, ds_times_conddp)
            dDdp = -dphidp + dphidp .* ( Status.Properties.S_1.*Status.Properties.cond_1 + ...
                Status.Properties.S_2.*Status.Properties.cond_2 ) + Status.Properties('Phi') .* ...
                ( ds_times_conddp(:,1) + ds_times_conddp(:,2) );
        end
        function dDdh = ComputeDDDh(obj, Status, ds_times_conddh)
            dDdh = Status.Properties('Phi') .* ( ds_times_conddh(:,1) + ds_times_conddh(:,2) );
        end

        % Other derivatives
        function drho_over_mudp = ComputeDrho_over_muDp(obj)
            drho_over_mudp = zeros(length(obj.Pindex),obj.NofPhases);
            for i=1:NofPhases
                rho_over_muTable = obj.TablePH.(['rho_over_mu_',num2str(i)]);
                drho_over_mudp(:,i) = obj.Phases(i).ComputeDrho_over_muDp(obj.Pindex, obj.Hindex, rho_over_muTable);
            end
        end
        function dUfdp = ComputeDUfDp(obj)
            UfTable = obj.TablePH.Uf;
            dUfdp = obj.Phases(1).ComputeDUfDp(obj.Pindex, obj.Hindex, UfTable);
        end
        function dUfdh = ComputeDUfDh(obj)
            UfTable = obj.TablePH.Uf;
            dUfdh = obj.Phases(1).ComputeDUfDh(obj.Pindex, obj.Hindex, UfTable);
        end        
        function dhdp = ComputeDhDp(obj)
            dhdp = zeros(length(obj.Pindex),obj.NofPhases);
            for i=1:NofPhases
                hTable = obj.TablePH.(['H_',num2str(i)]);
                dhdp(:,i) = obj.Phases(i).ComputeDhDp(obj.Pindex, obj.Hindex, hTable);
            end
        end
        function GetTotalFluidInternalEnergy(obj, Status)
            Uf = Status.Properties('Uf'); %%%%%%%%%%%%%%%%%%%%%%%%%
            Uf.Value = obj.TablePH.Uf( sub2ind(size(obj.TablePH.Uf), obj.Pindex, obj.Hindex) );
        end
        
 
    end
end
