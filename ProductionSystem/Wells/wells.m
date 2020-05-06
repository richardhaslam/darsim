% Wells 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef wells < handle
    properties
        Inj
        Prod
        NofInj
        NofProd
    end
    methods
        function AddInjector(obj, Injector)
            obj.Inj = [obj.Inj, Injector];
        end
        function AddProducer(obj, Producer)
            obj.Prod = [obj.Prod, Producer];
        end
        function InitializeFluxes(obj, n_phases, n_comp)
            % Injectors
            for i=1:obj.NofInj
                obj.Inj(i).QPhases = zeros(length(obj.Inj(i).Cells), n_phases);
                obj.Inj(i).Qh = zeros(length(obj.Inj(i).Cells), n_phases);
                obj.Inj(i).QComponents = zeros(length(obj.Inj(i).Cells), n_comp);
            end
            % Producers
            for i=1:obj.NofProd
                obj.Prod(i).QPhases = zeros(length(obj.Prod(i).Cells), n_phases);
                obj.Prod(i).Qh = zeros(length(obj.Prod(i).Cells), n_phases);
                obj.Prod(i).QComponents = zeros(length(obj.Prod(i).Cells), n_comp);
            end
        end
        function UpdateState(obj, Reservoir, FluidModel)
            K = Reservoir.K(:,1);
            switch FluidModel.name 
                case {'Geothermal_SinglePhase'} 
                    Mob = FluidModel.ComputePhaseMobilities(Reservoir.State.Properties('mu_1').Value);
                case {'Geothermal_MultiPhase'}
                    Mob = FluidModel.ComputePhaseMobilities(Reservoir.State);
                otherwise
                    Mob = FluidModel.ComputePhaseMobilities(Reservoir.State.Properties('S_1').Value);
            end
            
            % Injectors
            for i=1:obj.NofInj
                obj.Inj(i).UpdateState(Reservoir.State, K, FluidModel);         
            end
            % Producers
            for i=1:obj.NofProd
                obj.Prod(i).UpdateState(Reservoir.State, K, Mob, FluidModel);        
            end
        end
        function q = TotalFluxes(obj, Reservoir, Mobt)
            K = Reservoir.K(:,1);
            q = zeros(length(K), 1);
            
            % Injectors
            for i=1:obj.NofInj
                q = obj.Inj(i).TotalFlux(q, Reservoir.State.p, K);         
            end
            % Producers
            for i=1:obj.NofProd
                q = obj.Prod(i).TotalFlux(q, Reservoir.State.p, K, Mobt);        
            end
        end
    end
end