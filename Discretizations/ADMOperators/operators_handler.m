%  ADM operators handler base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 16 August 2016
%Last modified: 7 August 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef operators_handler < handle
    properties
        ProlongationBuilders
        ADMRest
        ADMProl
        FullOperatorsAssembler
    end
    methods
        function obj = operators_handler(cf)
            obj.ProlongationBuilders = prolongation_builder.empty;
        end
        function AddProlongationBuilder(obj, prolongationbuilder, index)
            obj.ProlongationBuilders(index) = prolongationbuilder;
        end
        function UpdateProlongationOperators(obj, ReservoirGrid, CoarseGrid, ProductionSystem)
            % Update basis functions for all variables
            for i=1:length(obj.ProlongationBuilders)
                obj.ProlongationBuilders(i).UpdateProlongationOperator(ReservoirGrid, CoarseGrid, ProductionSystem);
            end
        end
        function BuildADMOperators(obj, FineGrid, CoarseGrid, ADMGrid)
            start1 = tic; 
            % Restriction
            obj.ADMRestriction(ADMGrid, FineGrid);
            restriction = toc(start1);
            disp(['Restriction built in: ', num2str(restriction), ' s']);
            % Prolongation (do the pressure last coz it modifies ADMGrid)
            start2 = tic;
            for i=length(obj.ProlongationBuilders):-1:1
                obj.ADMProl{i} = obj.ProlongationBuilders(i).ADMProlongation(ADMGrid, FineGrid, CoarseGrid, obj.ADMRest);
            end
            prolongation = toc(start2);
            disp(['Prolongation built in: ', num2str(prolongation), ' s']);
        end
        function ADMRestriction(obj, ADMGrid, FineGrid)
              % Assemble dynamic FV restriction operator
              % Fine-scale cells
              obj.ADMRest = sparse(ADMGrid.Ntot, FineGrid.N);
              rows = 1:ADMGrid.N(1);
              columns = ADMGrid.CellIndex(1:ADMGrid.N(1))';
              obj.ADMRest(sub2ind(size(obj.ADMRest), rows, columns)) = 1;             
              % Coarse levels cells
              for c = ADMGrid.N(1) + 1:ADMGrid.Ntot 
                indexes = ADMGrid.GrandChildren{c};
                obj.ADMRest(c, indexes) = 1;
              end
        end
        function [Rest, Prol] = AssembleFullOperators(obj)
            [Rest, Prol] = obj.FullOperatorsAssembler.Assemble(obj.ADMRest, obj.ADMProl);
        end
    end
end