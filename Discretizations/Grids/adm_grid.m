%  dynamic grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 12 July 2016
%Last modified: 22 August 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef adm_grid < grid_darsim & matlab.mixin.Copyable
    properties
        Ntot
        level
        CellIndex
        MaxLevel
    end
    methods
        function Initialize(obj, n_total, numofactive, maxlevel)
            obj.MaxLevel = maxlevel;
            obj.N = numofactive;
            obj.Ntot = 0;
            obj.level = zeros(n_total, 1);
            obj.CoarseFactor = zeros(n_total, 3);
            obj.CellIndex = zeros(n_total, 1);
            obj.Fathers = zeros(n_total, maxlevel);
            obj.Children = cell(n_total, maxlevel);
            obj.Verteces = zeros(n_total, maxlevel);
        end
        function Update(obj, Nx, Nf, Nc, FineGrid)
            %       Nf     Nc
            %    --------------- 
            %    |      |      |
            % Nf |  I   |  0   | 
            %    |______|______|           
            %    |      |      |
            % Nx |      |      |
            %    |      |      |
            %    ---------------
            obj.N(:, obj.MaxLevel + 1) = 0;
            obj.MaxLevel = obj.MaxLevel - 1;
            obj.N(:, obj.MaxLevel + 1) = obj.N(:, obj.MaxLevel + 1) + Nx;
            obj.Ntot = sum(Nf) + sum(Nx);
            
            obj.CellIndex = [obj.CellIndex(1:sum(Nf)); zeros(sum(Nx), 1)];
            [~, MaxLevels] = size(obj.Fathers); 
            obj.Fathers = [obj.Fathers(1:sum(Nf),:); zeros(sum(Nx), MaxLevels)];
            obj.Verteces = [obj.Verteces(1:sum(Nf),:); zeros(sum(Nx), MaxLevels)];
            obj.level = [obj.level(1:sum(Nf)); ones(sum(Nx), 1)*obj.MaxLevel];
            
            % Make a new vector to store the new children
            ChildrenOfNc = obj.Children;
            obj.Children = vertcat( obj.Children(1:sum(Nf),:) , cell(sum(Nx), MaxLevels) );
            
            % Add the new cells to the new ADM grid
            h = sum(Nf) + 1;
            for CoarseNode = (sum(Nf)+1) : (sum(Nf)+sum(Nc))
                FineNodes = ChildrenOfNc{CoarseNode,1}; % if level l children belong to l-1
                n_children = length(FineNodes);
                obj.CellIndex(h:h + n_children-1) = FineNodes;
                for i = 1 : FineGrid.CoarseLevel
                    obj.Children(h:h + n_children-1, i) = FineGrid.Children(FineNodes, i);
                end
                obj.Fathers(h:h + n_children-1 , FineGrid.CoarseLevel+1:end ) = FineGrid.Fathers(FineNodes,:);
                obj.Verteces(h:h + n_children-1, FineGrid.CoarseLevel+1:end ) = FineGrid.Verteces(FineNodes,:);
                h = h + n_children;
            end
        end
    end
end