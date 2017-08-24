%  dynamic grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 12 July 2016
%Last modified: 22 August 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef adm_grid < grid_darsim
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
            obj.Children = cell(n_total, 1);
            obj.GrandChildren = cell(n_total, 1);
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
            obj.Fathers = [obj.Fathers(1:Nf,:); zeros(sum(Nx), MaxLevels)];
            obj.Verteces = [obj.Verteces(1:Nf,:); zeros(sum(Nx), MaxLevels)];
            obj.level = [obj.level(1:Nf); ones(sum(Nx), 1)*obj.MaxLevel];
            
            % Add the new cells
            h = sum(Nf) + 1;
            for CoarseNode = (sum(Nf)+1) : (sum(Nf)+sum(Nc))
                FineNodes = obj.Children{CoarseNode};
                n_children = length(FineNodes);
                obj.CellIndex(h:h + n_children-1) = FineNodes;
                obj.Fathers(h:h + n_children-1,:) = FineGrid.Fathers(FineNodes,:);
                obj.Verteces(h:h + n_children-1,:) = FineGrid.Verteces(FineNodes,:);
                h = h + n_children;
            end
        end
    end
end