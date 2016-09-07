%  ADM operators handler base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 16 August 2016
%Last modified: 16 August 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef operators_handler < handle
    properties
        R
        Pp
        ADMmap
        ADMRest
        ADMProlp
        ADMProls
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
             % Restriction
            obj.ADMRestriction(ADMGrid, FineGrid);
            
            % Prolongation
            obj.ADMProlongation(ADMGrid, FineGrid, CoarseGrid);
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
              obj.ADMRest = zeros(ADMGrid.Ntot, FineGrid.N);
              rows = 1:ADMGrid.N(1);
              columns = ADMGrid.CellIndex(1:ADMGrid.N(1))';
              obj.ADMRest(sub2ind(size(obj.ADMRest), rows, columns)) = 1;
              % Coarse levels cells
              for c = ADMGrid.N(1) + 1:ADMGrid.Ntot 
                obj.ADMRest(c, ADMGrid.GrandChildren{c}) = 1;
              end
              obj.ADMRest = sparse(obj.ADMRest);
        end
    end
    methods (Abstract)
        obj = BuildStaticOperators(obj);
    end
end
