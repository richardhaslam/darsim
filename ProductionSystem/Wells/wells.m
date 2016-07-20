% Wells 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 13 July 2016
%Last modified: 13 July 2016
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
                obj.Inj(i).QComponents = zeros(length(obj.Inj(i).Cells), n_comp);
            end
            % Producers
            for i=1:obj.NofProd
                obj.Prod(i).QPhases = zeros(length(obj.Prod(i).Cells), n_phases);
                obj.Prod(i).QComponents = zeros(length(obj.Prod(i).Cells), n_comp);
            end
        end
        function UpdateState(obj, Reservoir, FluidModel)
            K = Reservoir.K(:,1);
            Mob = FluidModel.ComputePhaseMobilities(Reservoir.State.S);
            
            % Injectors
            for i=1:obj.NofInj
                obj.Inj(i).UpdateState(Reservoir.State, K, FluidModel.NofPhases, FluidModel.NofComp);         
            end
            % Producers
            for i=1:obj.NofProd
                obj.Prod(i).UpdateState(Reservoir.State, K, Mob, FluidModel.NofPhases, FluidModel.NofComp);        
            end
        end
    end
end