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
        Ptable = (1:0.1:220)'; % The range of pressure for all the tables
        Htable = (20:1:4800)'; % The range of enthalpy for all the tables
        Pindex % The index of pressure value for lookup in the tables
        Hindex % The index of enthalpy value for lookup in the tables
        Pstepsize = 0.1;
        Hstepsize = 1; % implement this
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
            if obj.Pstepsize <= 1
                [~,obj.Pindex] = ismember(round(Status.Properties('P').Value,abs(log10(obj.Pstepsize))),round(obj.Ptable,abs(log10(obj.Pstepsize))));
            else 
                [~,obj.Pindex] = ismember(round(Status.Properties('P').Value,obj.Pstepsize),round(obj.Ptable,obj.Pstepsize));
            end
            if obj.Hstepsize <= 1
                [~,obj.Hindex] = ismember(round(Status.Properties('H').Value,abs(log10(obj.Hstepsize))),round(obj.Htable,abs(log10(obj.Hstepsize))));
            else
                [~,obj.Hindex] = ismember(round(Status.Properties('H').Value,obj.Hstepsize),round(obj.Htable,obj.Hstepsize));
            end               
%             [~,obj.Pindex] = ismember(round(Status.Properties('P').Value,1),round(obj.Ptable,1));
%             [~,obj.Hindex] = ismember(round(Status.Properties('H').Value,0),round(obj.Htable,0));
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
        function AddPhaseConductivities(obj, Status) % !!!
            cond = Status.Properties('cond_1');
            cond.Value = obj.Phases(1).AddConductivity(Status.Properties('P_1').Value, Status.Properties('T').Value);
        end

        function ComputePhaseDensities(obj, Status) 
            for i=1:obj.NofPhases
                rho = Status.Properties(['rho_', num2str(i)]);
                rhoTable = obj.TablePH.(['rho_', num2str(i)]);
                rho.Value = obj.Phases(i).ComputeDensity(obj.Pindex, obj.Hindex, rhoTable);
            end
        end
        function ComputePhaseSaturations(obj, Status)
            for i=1:obj.NofPhases
                S = Status.Properties(['S_',num2str(i)]); % S is just a pointer to a memory location here, so it has no value
                STable = obj.TablePH.(['S_',num2str(i)]); % |-> This is why you need num2str() again here
                S.Value = obj.Phases(i).ComputeSaturation(obj.Pindex, obj.Hindex, STable); 
            end
        end
        function ComputePhaseViscosities(obj, Status)
            for i=1:NofPhases
                mu = Status.Properties(['mu_',num2str(i)]);
                muTable = obj.TablePH.(['mu_',num2str(i)]);
                mu.Value = obj.Phases(i).ComputeViscosity(obj.Pindex, obj.Hindex, muTable);
            end
        end
        function ComputePhaseInternalEnergies(obj, Status)
            for i=1:NofPhases
                U = Status.Properties(['U_',num2str(i)]);
                UTable = obj.TablePH.(['U_',num2str(i)]);
                U.Value = obj.Phases(i).ComputeInternalEnergy(obj.Pindex, obj.Hindex, UTable);
            end
        end
        function ComputeTemperature(obj, Status)
            T = Status.Properties('T');
            T.Value = obj.TablePH.Temperature(sub2ind(size(obj.TablePH.Temperature), obj.Pindex, obj.Hindex));
        end
        function ComputeTotalDensity(obj, Status)   % In the accumulation term?
            % Compute the total density
            rhoT = Status.Properties('rhoT');
            rhoT.Value = Status.Properties('rho_1').Value .* Status.Properties('S_1').Value + ...
                Status.Properties('rho_2').Value .* Status.Properties('S_2').Value;
            % This can also be done with a table for the mixture density;
            % make sure both approaches are equal
        end

        function [drhodp, d2rhod2p] = ComputeDrhoDp(obj)
            drhodp = zeros(length(obj.Pindex),obj.NofPhases);
            d2rhod2p = zeros(length(obj.Pindex,obj.NofPhases));
            for i=1:obj.NofPhases
                rhoTable = obj.TablePH.(['rho_', num2str(i)]);
                [drhodp(:,i), d2rhod2p(:,i)] = obj.Phases(i).ComputeDrhoDp(obj.Pindex, obj.Hindex, rhoTable);
            end
        end   
        function [drhodh, d2rhod2h] = ComputeDrhoDh(obj)
            drhodh = zeros(length(obj.Hindex),obj.NofPhases);
            d2rhod2h = zeros(length(obj.Hindex,obj.NofPhases));
            for i=1:obj.NofPhases
                rhoTable = obj.TablePH.(['rho_', num2str(i)]);
                [drhodh(:,i), d2rhod2h(:,i)] = obj.Phases(i).ComputeDrhoDh(obj.Pindex, obj.Hindex, rhoTable);
            end
        end
        function [dSdp, d2Sd2p] = ComputeDSDp(obj)
            dSdp = zeros(length(obj.Pindex),obj.NofPhases);
            d2Sd2p = zeros(length(obj.Pindex,obj.NofPhases));
            for i=1:obj.NofPhases
                STable = obj.TablePH.(['S_', num2str(i)]);
                [dSdp(:,i), d2Sd2p(:,i)] = obj.Phases(i).ComputeDSDp(obj.Pindex, obj.Hindex, STable);
            end

        end
        function [dSdh, d2Sd2h] = ComputeDSDh(obj)
            dSdh = zeros(length(obj.Hindex),obj.NofPhases);
            d2Sd2h = zeros(length(obj.Hindex,obj.NofPhases));
            for i=1:obj.NofPhases
                STable = obj.TablePH.(['S_', num2str(i)]);
                [dSdh(:,i), d2Sd2h(:,i)] = obj.Phases(i).ComputeDSDh(obj.Pindex, obj.Hindex, STable);
            end
        end
        function [dmudp, d2mud2p] = ComputeDmuDp(obj)
            dmudp = zeros(length(obj.Pindex),obj.NofPhases);
            d2mud2p = zeros(length(obj.Pindex,obj.NofPhases));
            for i=1:NofPhases
                muTable = obj.TablePH.(['mu_',num2str(i)]);
                [dmudp(:,i),d2mud2p(:,i)] = obj.Phases(i).ComputeDmuDp(obj.Pindex, obj.Hindex, muTable);
            end
        end            
        function [dmudh, d2mud2h] = ComputeDmuDh(obj)
            dmudh = zeros(length(obj.Hindex),obj.NofPhases);
            d2mud2h = zeros(length(obj.Hindex,obj.NofPhases));
            for i=1:NofPhases
                muTable = obj.TablePH.(['mu_',num2str(i)]);
                [dmudh(:,i),d2mud2h(:,i)] = obj.Phases(i).ComputeDmuDh(obj.Pindex, obj.Hindex, muTable);
            end
        end
        function [dUdp, d2Ud2p] = ComputeDUDp(obj)
            dUdp = zeros(length(obj.Pindex),obj.NofPhases);
            d2Ud2p = zeros(length(obj.Pindex,obj.NofPhases));
            for i=1:NofPhases
                UTable = obj.TablePH.(['U_',num2str(i)]);
                [dUdp(:,i),d2Ud2p(:,i)] = obj.Phases(i).ComputeDUDp(obj.Pindex, obj.Hindex, UTable);
            end            
        end
        function [dUdh, d2Ud2h] = ComputeDUDh(obj)
            dUdh = zeros(length(obj.Hindex),obj.NofPhases);
            d2Ud2h = zeros(length(obj.Hindex,obj.NofPhases));
            for i=1:NofPhases
                UTable = obj.TablePH.(['U_',num2str(i)]);
                [dUdh(:,i),d2Ud2h(:,i)] = obj.Phases(i).ComputeDUDh(obj.Pindex, obj.Hindex, UTable);
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
               
        function Mob = ComputePhaseMobilities(obj, mu)
            Mob = zeros(length(obj.NofPhases),1);
            for i=1:NofPhases
                Mob(i,1) = 1./mu.Value;
            end
        end
        
        function dMobdp = ComputeDMobDp(obj,status)
            % For now, mobility has no dependency on pressure.
        end
        function v = ComputeVelocity(obj, Reservoir, mu)
%             virtual call
        end
    end
end