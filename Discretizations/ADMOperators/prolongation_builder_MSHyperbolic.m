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
        Pdelta
        Pdeltac
    end
    methods
        function obj = prolongation_builder_MSHyperbolic(n)
            obj@prolongation_builder(n)
            obj.P = cell(1, n);
            obj.Pdeltac = cell(1, n);
        end
        function BuildStaticOperators(obj, ProductionSystem, FluidModel, FineGrid, CrossConnections, maxLevel, CoarseGrid)
            % Build Restriction and Prolongation operators for static grids
            obj.R{1} = obj.MsRestriction(FineGrid, CoarseGrid(:, 1));
            obj.P{1} = obj.R{1}';
            for x = 2:maxLevel
                R = obj.MsRestriction(CoarseGrid(:, x-1), CoarseGrid(:,x));
                obj.R{x} = R*obj.R{x-1};
                obj.P{x} = obj.R{x}';
             end
        end
        function UpdateProlongationOperator(obj, FineGrid, CoarseGrid, ProductionSystem)
            %% Update the basis functions level by level (Update e = dSf/dSc)
            sf = ProductionSystem.Reservoir.State.Properties(obj.key).Value;
            sf0   =  ProductionSystem.Reservoir.State_old.Properties(obj.key).Value;
            for i=1:length(obj.R)
                obj.UpdateBasisFunctions(sf, sf0, i, CoarseGrid(i));
            end
        end
        function UpdateBasisFunctions(obj, sf, sf0, l, CoarseGrid)
            %% Update sat prolongation for level l
            % 1. Compute dSc from Sc^n = 1/V * sum(vSf^n)
            Sc  = obj.R{l}' * ((obj.R{l} * sf)  ./ sum(obj.R{l}, 2));
            Sc0 = obj.R{l}' * ((obj.R{l} * sf0) ./ sum(obj.R{l}, 2));
            deltac = Sc - Sc0;
            deltaf = sf - sf0;
            epsilon = deltaf ./ deltac;
            epsilon(isnan(epsilon)) = 1;
            epsilon(epsilon == Inf) = 1;
            epsilon(abs(deltac) < 1e-4) = 1;
            epsilon(abs(deltaf) < 1e-4) = 1;
            % 2. Update the prolongation operator
            obj.Pdelta  = deltaf;
            obj.Pdeltac{l} = deltac;
            for c=1:CoarseGrid.N
                fsI = CoarseGrid.GrandChildren(c,:);
                %obj.P{l}(fsI, c) = epsilon(fsI);
                obj.P{l}(fsI, c) = epsilon(fsI)/max(epsilon(fsI));
            end
        end
        function ADMProl = ADMProlongation(obj, ADMGrid, FineGrid, CoarseGrid, ADMRest)
            ADMProl = ADMRest';
            % Coarse levels cells
            for c = ADMGrid.N(1) + 1:ADMGrid.Ntot
                indexes = ADMGrid.GrandChildren{c};
                ADMProl(indexes, c) = obj.P{ADMGrid.level(c)}(indexes, ADMGrid.CellIndex(c));
            end
        end
    end
end