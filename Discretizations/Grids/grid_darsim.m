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
        Children
        GrandChildren
        Verteces
        CoarseFactor
        GridCoords
        ActiveTime
        
    end
    methods
        function Initialise(obj, maxLevel)
            % They can be very heterogeneous so I use cell arrays
            obj.Fathers = zeros(obj.N, maxLevel);
            obj.Children = cell(obj.N, 1);
            obj.GrandChildren = cell(obj.N, 1);
            obj.Verteces = zeros(obj.N, maxLevel);
            obj.CoarseFactor = zeros(obj.N, 3);            
        end
        function CopyGridEntries(obj, Grid, Nc_global, level)
            for c = 1:Grid.N
                obj.Fathers(c + Nc_global(level), :) = Grid.Fathers(c, :) + Nc_global(2:end);
                obj.Children{c + Nc_global(level)} = Grid.Children(c,:) + Nc_global(level);
                obj.GrandChildren{c + Nc_global(level)} = Grid.GrandChildren(c,:) + Nc_global(1);
                obj.Verteces(c + Nc_global(level), :) = Grid.Verteces(c,:);
                obj.CoarseFactor(c + Nc_global(level), :) = Grid.CoarseFactor;
            end
        end
    end
end