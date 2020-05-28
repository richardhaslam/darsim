%  grid base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 12 July 2016
%Last modified: 22 August 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef grid_darsim < matlab.mixin.Heterogeneous & handle
    properties
        N
        Depth
        Active
        Fathers
        Children = {}
        Verteces
        CoarseFactor
        GridCoords
        ActiveTime
        %CoarseNodeIndex
        CoarseLevel
    end
    methods
        function Initialise(obj, currentLevel, maxLevel)
            % They can be very heterogeneous so I use cell arrays
            obj.Children = cell(obj.N, currentLevel);
            obj.Fathers = zeros(obj.N, maxLevel-currentLevel);
            obj.Verteces = zeros(obj.N, maxLevel-currentLevel);
            obj.CoarseFactor = zeros(obj.N, 3);
        end
        function CopyGridEntries(obj, Grid, Nc_global, level)
            obj.Fathers = [obj.Fathers; Grid.Fathers];
            obj.Children
            for c = 1:Grid.N
                if ~isempty(Grid.Fathers)
                    obj.Fathers(c + Nc_global(2+level), :) = Grid.Fathers(c, :) + Nc_global(2+level:end);
                end
                obj.Children{c + Nc_global(level)} = Grid.Children{c,:} + Nc_global(level);
                if ~isempty(Grid.Verteces)
                    obj.Verteces(c + Nc_global(level), :) = Grid.Verteces(c,:);
                    obj.CoarseFactor(c + Nc_global(level), :) = Grid.CoarseFactor;
                end
            end
        end
    end
end