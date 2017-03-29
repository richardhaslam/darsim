%  ADM operators handler base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 16 August 2016
%Last modified: 27 September 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef operators_handler < handle
    properties
        R
        Pp
        ADMmap
        ADMRest
        ADMProlp
        ADMProls
        ADMProlAv
        FullOperatorsAssembler
    end
    methods
        function obj = operators_handler(n, CF)
            obj.R = cell(1, n);
            obj.Pp = cell(1, n);
            obj.ADMmap = adm_map(CF);
        end
        function MsR = MsRestriction(obj, FineGrid, CoarseGrid)
            Nf = FineGrid.N;
            Nc = CoarseGrid.N;
            %% MSFV Restriction Operator
            MsR = zeros(Nc, Nf);
            for c = 1:Nc
                MsR(c, CoarseGrid.Children(c,:)) = 1;
            end
            MsR = sparse(MsR);
        end
        function BuildADMOperators(obj, FineGrid, CoarseGrid, ADMGrid)
            start1 = tic; 
            % Restriction
            obj.ADMRestriction(ADMGrid, FineGrid);
            restriction = toc(start1);
            disp(['Restriction built in: ', num2str(restriction), ' s']);
            % Prolongation
            start2 = tic;
            obj.ADMProlongation(ADMGrid, FineGrid, CoarseGrid);
            prolongation = toc(start2);
            disp(['Prolongation built in: ', num2str(prolongation), ' s']);
        end
        function ADMRestriction(obj, ADMGrid, FineGrid)
%             Rf = speye(obj.ADMmap.Nf);
%             Rc = zeros(obj.ADMmap.Nc, obj.ADMmap.Nx);
%             % Restriction operator
%             for c = 1:obj.ADMmap.Nc
%                     % Active coarse cells
%                     i = (obj.ADMmap.Nf+1) * c;
%                     f = obj.ADMmap.Nf + CF;
%                     Rc(c, i:f) = 1;
%             end
%             Rest = [Rf, zeros(obj.ADMmap.Nf, obj.ADMmap.Nx);
%                 zeros(obj.ADMmap.Nc, obj.ADMmap.Nx), Rc];
%             obj.ADMRest{level} = sparse(Rest);
              
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              % Fine-scale cells
              
              obj.ADMRest = sparse(ADMGrid.Ntot, FineGrid.N);
              rows = 1:ADMGrid.N(1);
              columns = ADMGrid.CellIndex(1:ADMGrid.N(1))';
              obj.ADMRest(sub2ind(size(obj.ADMRest), rows, columns)) = 1;
              obj.ADMProlAv = obj.ADMRest';              
              % Coarse levels cells
              for c = ADMGrid.N(1) + 1:ADMGrid.Ntot 
                indexes = ADMGrid.GrandChildren{c};
                obj.ADMRest(c, indexes) = 1;
                obj.ADMProlAv(indexes, c) = 1/prod(ADMGrid.CoarseFactor(c,:));
              end
        end
        function [Rest, Prol] = AssembleFullOperators(obj)
            [Rest, Prol] = obj.FullOperatorsAssembler.Assemble(obj.ADMRest, obj.ADMProlp, obj.ADMProls);
        end
    end
    methods (Abstract)
        obj = BuildStaticOperators(obj);
    end
end
