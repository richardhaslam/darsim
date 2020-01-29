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
        function GetTableIndex(obj, p, h)
            [~,obj.Pindex] = ismember(round(p,1),round(obj.Ptable,1));
            [~,obj.Hindex] = ismember(round(h,0),round(obj.Htable,0));
        end
        function InitializeInjectors(obj, Inj) % This one needs to change to something that looks like Init_Inj in SinglePhase
            for i=1:length(Inj)
                Inj(i).z = [1 0];
                Inj(i).x = [1 0 0 1];
                Inj(i).S = 1;
                for ph=1:obj.NofPhases
                    Inj(i).rho(:, ph)= obj.Phases(ph).ComputeDensity(Inj(i).p);
                end
                Inj(i).x2 = 1 - Inj(i).x;
                Inj(i).Mob = obj.ComputePhaseMobilities(Inj(i).S);   
            end
        end   % is this one correct?
        
%         function ComputeTemperature()
%         end
        
        function ComputePhaseDensities(obj, Status) % correct
            for i=1:obj.NofPhases
                rho = Status.Properties(['rho_', num2str(i)]);
                rhoTable = TablePH.(['rho_', num2str(i)]);
                rho.Value = obj.Phases(i).ComputeDensity(obj.Pindex, obj.Hindex, rhoTable);
            end
        end
        function ComputePhaseSaturations(obj, Status)
            for i=1:obj.NofPhases
                S = Status.Properties(['S_',num2str(i)]);
                TemporaryPhaseIndex = TablePH.(S{1});
                S.Value = TemporaryPhaseIndex(sub2ind(size(TemporaryPhaseIndex), obj.Pindex, obj.Hindex));
            end
        end
        function AddPhaseConductivities(obj, Status) % !!!
            cond = Status.Properties('cond_1');
            cond.Value = obj.Phases(1).AddConductivity(Status.Properties('P_1').Value, Status.Properties('T').Value);
        end
        function ComputePhaseInternalEnergies(obj, Status)
            for i=1:NofPhases
                U = Status.Properties(['U_1',num2str(i)]);
                TemporaryPhaseIndex = TablePH.(U{1});
                U.Value = TemporaryPhaseIndex(sub2ind(size(TemporaryPhaseIndex), obj.Pindex, obj.Hindex));
            end
        end
        function ComputePhaseViscosities(obj, Status)
            for i=1:NofPhases
                mu = Status.Properties(['mu_1',num2str(i)]);
                TemporaryPhaseIndex = TablePH.(mu{1});
                mu.Value = TemporaryPhaseIndex(sub2ind(size(TemporaryPhaseIndex), obj.Pindex, obj.Hindex));
            end
        end


        function [drhodp, d2rhod2p] = ComputeDrhoDp(obj)
            drhodp = zeros(length(obj.Pindex,obj.NofPhases));
            d2rhod2p = zeros(length(obj.Pindex,obj.NofPhases));
            for i=1:obj.NofPhases
                rhoTable = TablePH.(['rho_', num2str(i)]);
                [drhodp(:,i), d2rhod2p(:,i)] = obj.Phases(i).ComputeDrhoDp(obj.Pindex, obj.Hindex, rhoTable);
            end
        end
        
        function [drhodh, d2rhod2h] = ComputeDrhoDh(obj)
            drhodh = zeros(length(obj.Hindex,obj.NofPhases));
            d2rhod2h = zeros(length(obj.Hindex,obj.NofPhases));
            for i=1:obj.NofPhases
                rhoTable = TablePH.(['rho_', num2str(i)]);
                [drhodh(:,i), d2rhod2h(:,i)] = obj.Phases(i).ComputeDrhoDh(obj.Pindex, obj.Hindex, rhoTable);
            end
        end

%         function ComputeTotalDensity(obj, Status)   % In the accumulation term?
%             % Compute the total density
%             rhoT = Status.Properties('rhoT');
%             rhoT.Value = Status.Properties('rho_1').Value; 
%             % For 1 phase rhoT is rho1
%         end
               
        
        function [dmudT,d2mudT2] = ComputeDmuDT(obj, Status)
                [~,dmudT,d2mudT2] = obj.Phases.ComputeViscosity(Status.Properties('T').Value);
        end
        function Mob = ComputePhaseMobilities(obj, mu) 
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
        function v = ComputeVelocity(obj, Reservoir, mu)
%             virtual call
        end
    end
end