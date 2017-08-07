%  Prolongation builder MS for hyperbolic variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 4 August 2017
%Last modified: 7 August 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef prolongation_builder_MSHyperbolic < prolongation_builder
    properties
        key = 'S_1';
    end
    methods
        function obj = prolongation_builder_MSHyperbolic(n, cf)
            obj@prolongation_builder(n)
            obj.P = cell(1, n);
        end
        function BuildStaticOperators(obj, CoarseGrid, FineGrid, maxLevel, K, s, FluidModel)
            % Build Restriction and Prolongation operators for static grids
            obj.R{1} = MsRestriction(obj, FineGrid, CoarseGrid(1));
            obj.P{1} = obj.R{1}';
            for x = 2:maxLevel
                R = MsRestriction(obj, CoarseGrid(x-1), CoarseGrid(x));
                obj.R{x} = R*obj.R{x-1};
                obj.P{x} = obj.R{x}';
             end
        end
        function UpdateProlongationOperator(obj, FineGrid, CoarseGrid, ProductionSystem)
            %% Update the basis functions level by level (Update e = dSf/dSc)
            sf = ProductionSystem.Reservoir.State.Properties(obj.key).Value;
            sf0   =  ProductionSystem.Reservoir.State_old.Properties(obj.key).Value;
            for i=1:length(obj.R)
                obj.UpdateBasisFunctions(sf, sf0, i);
            end
        end
        function UpdateBasisFunctions(obj, sf, sf0, l)
            %% Update sat prolongation for level l
            % 1. Compute dSc from Sc^n = 1/V * sum(vSf^n)
            Sc  = obj.R{l}' * ((obj.R{l} * sf)  ./ sum(obj.R{l}, 2));
            Sc0 = obj.R{l}' * ((obj.R{l} * sf0) ./ sum(obj.R{l}, 2));
            deltac = Sc - Sc0;
            deltaf = sf - sf0;
            epsilon = deltaf ./ deltac;
            % 2. Update the prolongation operator
            
            obj.P{l} = xx;
        end
        function ADMProl = ADMProlongation(obj, ADMGrid, FineGrid, CoarseGrid, ADMRest)
            ADMProl = ADMRest';
            % Coarse levels cells
            for c = ADMGrid.N(1) + 1:ADMGrid.Ntot
                indexes = ADMGrid.GrandChildren{c};
                ADMProl(indexes, c) = obj.P{ADMGrid.level{c}}(indexes, c);
            end
        end
    end
end