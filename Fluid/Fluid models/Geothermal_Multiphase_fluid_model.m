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
        SinglePhase % is this necessary ??
    end
    methods
        function obj = Geothermal_Multiphase_fluid_model(n_phases)
            obj@fluid_model(n_phases, n_phases);
            obj.name = 'Geothermal_Multiphase';
        end
        function SinglePhase = Flash(obj, Status)
            % Composition in this case is fixed to be 1 and 0
            SinglePhase = zeros(length(Status.Properties('S_1').Value), 1);
            SinglePhase (Status.Properties('S_1').Value == 1) = 1;
            SinglePhase (Status.Properties('S_2').Value == 1) = 2;            
        end
        % What does 'SinglePhase' do ??
        
        
        % do we need GetTableIndex function here as well??
        
        
        function InitializeInjectors(obj, Inj)
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
        
        function ComputePhaseDensities(obj, Status)
            for i=1:obj.NofPhases
                rho = Status.Properties(['rho_', num2str(i)]);
                TemporaryPhaseIndex = TablePH.(rho{1});
                rho.Value = TemporaryPhaseIndex(sub2ind(size(TemporaryPhaseIndex), obj.Pindex, obj.Hindex));
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


        function ComputeDensityDerivatives(obj, Status)
            for i=1:obj.NofPhases
                rho = Status.Properties(['rho_', num2str(i)]); % only for structure table_index in this case (line 'TemporaryPhaseIndex')
                
                drhodp = Status.Properties(['drhodp_',num2str(i)]);
                d2rhod2p = Status.Properties(['d2rhod2p_',num2str(i)]);
                TemporaryPhaseIndex = TablePH.(rho{1});
                % 1st derivative 
                [~,table_drhodp] = gradient(TemporaryPhaseIndex,1,0.1); %gradient('matrix','stepsize hor(j)','stepsize vert(i)') and you have P,H as i,j
                drhodp.Value = table_drhodp(sub2ind(size(table_drhodp), obj.Pindex, obj.Hindex));

                % 2nd derivative
                [~,table_d2rhod2p] = gradient(table_drhodp,0.1); % specify stepsize for pressure (make this generic)
                d2rhod2p.Value = table_d2rhod2p(sub2ind(size(table_d2rhod2p), obj.Pindex, obj.Hindex));
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