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
        function obj = operators_handler(n)
            obj.R = cell(1, n);
            obj.Pp = cell(1, n);
            obj.ADMmap = adm_map();
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
            obj.ADMRestriction(ADMGrid);
            
            %Other Levels
            for i=ADMGrid.MaxLevel:-1:2
                % Prolongation
                obj.ADMProlongation(i);
            end
            % FIRST LEVEL
            % Prolongations
            obj.ADMProlongation(1);
        end
        function ADMRestriction(obj, ADMGrid)
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
              obj.ADMRest(1:ADMGrid.N(1), ADMGrid.CellIndex(1:ADMGrid.N(1))) = 1;
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
