%  ADM map 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef adm_map < handle
    properties
       Nf % all cells belonging to all levels l < lc
       Nc % all cells belonging to level lc
       Nx % number of all the cells of level "lc - 1" generated from lc
       CF
       Verteces 
       OriginalIndexVerteces
       OriginalIndexNc
       OriginalIndexNx
       OriginalIndexNf
    end
    methods
        function obj = adm_map(cf)
            obj.CF = cf;
        end
        function Update(obj, ADMGrid, FineGrid, level)
            % ADMGrid.N contains the number of the grid cells at each level for ech medium for the current ADMGrid
            % Structure:
            %
            % Medium Index  ---  NumOfGrids at Level 0 (Finescale)  ---  NumOfGrids at Level 1  ---    NumOfGrids at Level 2  ---   NumOfGrids at Level 3  --- ...
            % Reservoir                    N_res_L0                           N_res_L1                      N_res_L2                     N_res_L3
            % Fracture 1                  N_frac1_L0                         N_frac1_L1                    N_frac1_L2                   N_frac1_L3
            % Fracture 2                  N_frac2_L0                         N_frac2_L1                    N_frac2_L2                   N_frac2_L3
            % Fracture 3                  N_frac3_L0                         N_frac3_L1                    N_frac3_L2                   N_frac3_L3
            %    .
            %    .
            %    .
            
            % The "FineGrid" is the global grid (reservoir and fractures combined) belonging to "level-1"
            
            % CF: the coarsening factor (x.y.z multiplied) for each medium
            obj.CF = obj.CF(1:size(ADMGrid.N,1));
            % Nf: Number of all the cells at each level lower than current level ("level") (reservoir and fractures combined)
            obj.Nf = sum(ADMGrid.N(:,1:level), 1);
            % NC: Number of all the cells at current level ("level") (separate for reservoir and each fracture)
            obj.Nc = ADMGrid.N(:, level+1);
            
            % NumOfChildrenFull: Number of children at each level belonging to each of the ADM grid cell
            NumOfChildrenFull = cellfun('length',ADMGrid.Children);
            % Looping over each medium (reservor and fractures)
            for m = 1 : size(ADMGrid.N,1)
                % Finding the start and end index of this medium from the ADMGrid
                Start = sum(sum(ADMGrid.N(1:m-1,:)))+1;
                End = sum(sum(ADMGrid.N(1:m,:)));
                % Extracting the NumOfChildren for this medium
                NumOfChildrenMedia = NumOfChildrenFull(Start:End, :);
                % Nx = number of finer grid cells inside the grid cells of the medium at current level ("level")
                % As the "Children" contains all the childeren at all levels, we use the first column to take
                % the children from only one level finer and not more.
                obj.Nx(m,1) = sum( NumOfChildrenMedia( ADMGrid.level(Start:End)==level , 1 ) );
            end
            
            % Original Indeces of Nf
            obj.OriginalIndexNf = ADMGrid.CellIndex(1:sum(obj.Nf));
            
            % Verteces
            % obj.Verteces = find(ADMGrid.Verteces(:,level)); from the old implementation. Replacing it with the line below.
            obj.Verteces = find(all(ADMGrid.Verteces(:,1:level)~=0, 2));
            obj.OriginalIndexVerteces = ADMGrid.Fathers(obj.Verteces, level);
            SelfGridsOfThisLevel = find(all(ADMGrid.Verteces(:,1:level)==-1, 2));
            obj.OriginalIndexVerteces(obj.OriginalIndexVerteces==-1) = ADMGrid.CellIndex(SelfGridsOfThisLevel);
            
            % Coarse nodes that will be prolonged
            obj.OriginalIndexNc = ADMGrid.CellIndex(sum(obj.Nf)+1 : end);
            
            % Update ADMGrid
            if level > 1
                ADMGrid.Update(obj.Nx, obj.Nf, obj.Nc, FineGrid);
                
                % New Fine-Grid indeces
                obj.OriginalIndexNx = ADMGrid.CellIndex(sum(obj.Nf)+1 : end);
            end
        end
    end
end