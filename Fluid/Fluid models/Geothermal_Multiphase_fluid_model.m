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
        Ptable = (1:0.1:220)'.*1e5; % The range of pressure for all the tables
        Htable = (20:1:4800)'*1e3; % The range of enthalpy for all the tables
        Pindex % The index of pressure value for lookup in the tables
        Hindex % The index of enthalpy value for lookup in the tables
        Pstepsize = 1e4;
        Hstepsize = 1e3; % implement this
    end
    methods
        function obj = Geothermal_Multiphase_fluid_model(n_phases)
            obj@fluid_model(n_phases, n_phases);
            obj.name = 'Geothermal_Multiphase';
        end
        function SinglePhase = Flash(obj, Status)
            % What does 'SinglePhase' do ??
            % Composition in this case is fixed to be 1 and 0
            SinglePhase = zeros(length(Status.Properties('S_1').Value), 1);
            SinglePhase (Status.Properties('S_1').Value == 1) = 1;
            SinglePhase (Status.Properties('S_2').Value == 1) = 2;            
        end
        function GetTableIndex(obj, Status)
            [~,obj.Pindex] = ismember( round(Status.Properties('P').Value,-log10(obj.Pstepsize)), round(obj.Ptable,-log10(obj.Pstepsize)) );
            [~,obj.Hindex] = ismember( round(Status.Properties('H').Value,-log10(obj.Hstepsize)), round(obj.Htable,-log10(obj.Hstepsize)) );
%             [~,obj.Pindex] = ismember(round(Status.Properties('P').Value,-4),round(obj.Ptable,-4));
%             [~,obj.Hindex] = ismember(round(Status.Properties('H').Value,-3),round(obj.Htable,-3));
        end

        function GetPhaseDensities(obj, Status) 
            for i=1:obj.NofPhases
                rho = Status.Properties(['rho_', num2str(i)]);
                rhoTable = obj.TablePH.(['rho_', num2str(i)]);
                rho.Value = obj.Phases(i).GetDensity(obj.Pindex, obj.Hindex, rhoTable);
            end
        end
        function GetPhaseSaturations(obj, Status)
            for i=1:obj.NofPhases
                S = Status.Properties(['S_',num2str(i)]); % S is just a pointer to a memory location here, so it has no value
                STable = obj.TablePH.(['S_',num2str(i)]); % |-> This is why you need num2str() again here
                S.Value = obj.Phases(i).GetSaturation(obj.Pindex, obj.Hindex, STable); 
            end
        end
        function GetPhaseViscosities(obj, Status)
            for i=1:NofPhases
                mu = Status.Properties(['mu_',num2str(i)]);
                muTable = obj.TablePH.(['mu_',num2str(i)]);
                mu.Value = obj.Phases(i).GetViscosity(obj.Pindex, obj.Hindex, muTable);
            end
        end
        function GetPhaseInternalEnergies(obj, Status)
            for i=1:NofPhases
                U = Status.Properties(['U_',num2str(i)]);
                UTable = obj.TablePH.(['U_',num2str(i)]);
                U.Value = obj.Phases(i).GetInternalEnergy(obj.Pindex, obj.Hindex, UTable);
            end
        end
        function GetPhaseThermalConductivities(obj, Status) 
            for i=1:NofPhases
                ThermCond = Status.Properties(['cond_',num2str(i)]);
                ThermCondTable = obj.TablePH.(['cond_',num2str(i)]);
                ThermCond.Value = obj.Phases(i).GetConductivity(obj.Pindex, obj.Hindex, ThermCondTable);
            end
        end
        function GetPhaseEnthalpies(obj, Status)
            for i=1:NofPhases
                PhaseEnthalpy = Status.Properties(['H_',num2str(i)]);
                PhaseEnthalpyTable = obj.TablePH.(['H_',num2str(i)]);
                PhaseEnthalpy.Value = obj.Phases(i).GetPhaseEnthalpy(obj.Pindex, obj.Hindex, PhaseEnthalpyTable);    
            end
        end
        
        function GetTemperature(obj, Status)
            T = Status.Properties('T');
            T.Value = obj.TablePH.Temperature(sub2ind(size(obj.TablePH.Temperature), obj.Pindex, obj.Hindex));
        end
        function GetTotalDensity(obj, Status)   
            rhoT = Status.Properties('rhoT');
            rhoT.Value = obj.TablePH.rhoT( sub2ind(size(obj.TablePH.rhoT), obj.Pindex, obj.Hindex) ); 
        end
        function GetTotalFluidInternalEnergy(obj, Status)
            Uf = Status.Properties('Uf'); %%%%%%%%%%%%%%%%%%%%%%%%%
            Uf.Value = obj.TablePH.Uf( sub2ind(size(obj.TablePH.Uf), obj.Pindex, obj.Hindex) );
        end
        
%         function ComputeThermalConductivityTensor(obj, Status)
%             D = Status.Properties('D');
%             D.Value = (1 - Status.Properties('Phi')) .* D_rock + Status.Properties('Phi') .* ...
%                 ( Status.Properties('cond_1') .* Status.Properties('S_1') + Status.Properties('cond_2') .* Status.Properties('S_2') );
%         end
        
        
        % Derivatives mass balance
        function drhodp = ComputeDrhoDp(obj)
            drhodp = zeros(length(obj.Pindex),obj.NofPhases);
            for i=1:obj.NofPhases
                rhoTable = obj.TablePH.(['rho_', num2str(i)]);
                drhodp(:,i) = obj.Phases(i).ComputeDrhoDp(obj.Pindex, obj.Hindex, rhoTable);
            end
        end
        function drhodh = ComputeDrhoDh(obj)
            drhodh = zeros(length(obj.Hindex),obj.NofPhases);
            for i=1:obj.NofPhases
                rhoTable = obj.TablePH.(['rho_', num2str(i)]);
                drhodh(:,i) = obj.Phases(i).ComputeDrhoDh(obj.Pindex, obj.Hindex, rhoTable);
            end
        end
        function drhoTdp = ComputeDrhoTDp(obj)
%             drhoTdp = zeros(length(obj.Pindex),1);
            rhoTTable = obj.TablePH.rhoT;
            drhoTdp = obj.Phases(1).ComputeDrhoTDp(obj.Pindex, obj.Hindex, rhoTTable);
        end
        function drhoTdh = ComputeDrhoTDh(obj)
%             drhoTdh = zeros(length(obj.Hindex),1);
            rhoTTable = obj.TablePH.rhoT;
            drhoTdh = obj.Phases(1).ComputeDrhoTDh(obj.Pindex, obj.Hindex, rhoTTable);
        end
        function drho_over_mudp = ComputeDrho_over_muDp(obj)
            drho_over_mudp = zeros(length(obj.Pindex),obj.NofPhases);
            for i=1:NofPhases
                rho_over_muTable = obj.TablePH.(['rho_over_mu_',num2str(i)]);
                drho_over_mudp(:,i) = obj.Phases(i).ComputeDrho_over_muDp(obj.Pindex, obj.Hindex, rho_over_muTable);
            end
        end

        % Derivatives energy balance
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
        function drho_times_hdp = ComputeDrho_times_hDp(obj)
            drho_times_hdp = zeros(length(obj.Pindex),obj.NofPhases);
            for i=1:NofPhases
                rho_times_hTable = obj.TablePH.(['rho_times_H_',num2str(i)]);
                drho_times_hdp(:,i) = obj.Phases(i).ComputeDrho_times_hDp(obj.Pindex, obj.Hindex, rho_times_hTable);
            end
        end        
        function drho_times_hdh = ComputeDrho_times_hDh(obj)
            drho_times_hdh = zeros(length(obj.Hindex),obj.NofPhases);
            for i=1:NofPhases
                rho_times_hTable = obj.TablePH.(['rho_times_H_',num2str(i)]);
                drho_times_hdh(:,i) = obj.Phases(i).ComputeDrho_times_hDh(obj.Pindex, obj.Hindex, rho_times_hTable);
            end
        end        
        function [dTdp, d2Td2p] = ComputeDTDp(obj)
            TTable = obj.TablePH.Temperature; % What if we just pass obj.TablePH.Temperature to the function ??
            [dTdp, d2Td2p] = obj.Phases(1).ComputeDTDp(obj.Pindex, obj.Hindex, TTable);
        end
        function [dTdh, d2Td2h] = ComputeDTDh(obj)
            TTable = obj.TablePH.Temperature;
            [dTdh, d2Td2h] = obj.Phases(1).ComputeDTDh(obj.Pindex, obj.Hindex, TTable);
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
        
        

        function Mob = ComputePhaseMobilities(obj, mu)
            Mob = zeros(length(obj.NofPhases),1);
            for i=1:NofPhases
                Mob(i,1) = 1./mu.Value;
            end
        end       
        function dMobdp = ComputeDMobDp(obj,status)
            % For now, mobility has no dependency on pressure.
        end
        
        function InitializeInjectors(obj, Inj)
            for i=1:length(Inj)
                Inj(i).z = 1;
                Inj(i).x = [1 0];
                Inj(i).S = 1;
                for ph=1:obj.NofPhases
                    rhoTable = obj.TablePH.(['rho_', num2str(i)]);
                    muTable = obj.TablePH.(['mu_1',num2str(i)]);
                    Inj(i).rho(:, ph)= obj.Phases(ph).ComputeDensity(Inj(i).p, Inj(i).h, rhoTable);
%                     Inj(i).h(:, ph)= obj.Phases(ph).ComputeEnthalpy(Inj(i).p, Inj(i).T); % if we initialize with temperature, this one is still necessary
                    mu = obj.Phases(ph).ComputeViscosity(Inj(i).p, Inj(i).h, muTable);
                end
                Inj(i).Mob = 1/mu;   
            end
        end
        
        function v = ComputeVelocity(obj, Reservoir, mu)
%             virtual call
        end          
    end
end
