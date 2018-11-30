% OBL Formulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Jeremy Yong
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef OBL_formulation < fim_formulation
    properties
        Table1
        Table2
        Table3
        Table4
    end
    methods
        function obj = OBL_formulation(n_components)
            obj@fim_formulation(n_components);
        end
        function Reset(obj)
            % It s a virtual call
        end
        function SavePhaseState(obj)
            % It s a virtual call
        end
        function CreateTables(obj, File)
            % Here you read tables
            obj.Table1 = % read from ADGPRS
            obj.Table2 = % read from ADGPRS
        end
        function ComputePropertiesAndDerivatives(obj, ProductionSystem, FluidModel)
            % Jeremy does interpolation here.
            % Check Immiscible_formulation.m 
            % you will need to access obj.Table1 etc... to create obj.Tph1 
            
        end
        function Residual = BuildResidual(obj,ProductionSystem, DiscretizationModel, dt, State0)
            % Jeremy builds Residual here
            
        end
        function Jacobian = BuildJacobian(obj, ProductionSystem, DiscretizationModel, dt)
            % Jeremy builds Jacobian!!
            
            % BUILD FIM JACOBIAN BLOCK BY BLOCK
           
           
        end
        function [Jwp, JwS, Jnwp, JnwS] = AddWellsToJacobian(obj, Jwp, JwS, Jnwp, JnwS, State, Wells, K)
            Inj = Wells.Inj;
            Prod = Wells.Prod;
            %Injectors
            for i=1:length(Inj)
                a = Inj(i).Cells;
                for j=1:length(a)
                    Jnwp(a(j),a(j)) = Jnwp(a(j),a(j)) + Inj(i).PI*K(a(j))*Inj(i).Mob(:, 2)*Inj(i).rho(:, 2);
                    Jwp(a(j),a(j)) = Jwp(a(j),a(j)) + Inj(i).PI*K(a(j))*Inj(i).Mob(:, 1)*Inj(i).rho(:, 1);
                end
            end
            %Producers
            for i=1:length(Prod)
                b = Prod(i).Cells;
                for j=1:length(b)
                    Jnwp(b(j),b(j)) = Jnwp(b(j),b(j)) + Prod(i).PI*K(b(j)).*obj.Mob(b(j), 2) .* State.rho(b(j), 2)...
                     - Prod(i).PI * K(b(j)) * obj.Mob(b(j), 2) * obj.drhodp(b(j),2) * (Prod(i).p - State.p(b(j)));                    
                    Jwp(b(j),b(j)) = Jwp(b(j),b(j)) + Prod(i).PI*K(b(j)).*obj.Mob(b(j), 1) .* State.rho(b(j), 1)...
                     - Prod(i).PI * K(b(j)) * obj.Mob(b(j), 1) * obj.drhodp(b(j),1) * (Prod(i).p - State.p(b(j)));
                    
                    JwS(b(j),b(j)) = JwS(b(j),b(j)) - Prod(i).PI*K(b(j)).* State.rho(b(j), 1) .* (Prod(i).p - State.p(b(j))).*obj.dMob(b(j), 1);
                    JnwS(b(j),b(j)) = JnwS(b(j),b(j)) - Prod(i).PI*K(b(j)).* State.rho(b(j), 2) .* (Prod(i).p - State.p(b(j))).*obj.dMob(b(j), 2);
                end
            end
        end
        function UpdateState(obj, delta, Status, FluidModel)
           
        end
    end
end