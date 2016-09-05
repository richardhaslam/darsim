%  dynamic grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 12 July 2016
%Last modified: 12 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef adm_grid < grid
    properties
        Ntot
        I
        J
        level
        CellIndex
        MaxLevel
    end
    methods
        function Initialize(obj, Ntotal, NumberOfActive, maxlevel)
            obj.MaxLevel = maxlevel;
            obj.N = NumberOfActive; 
            obj.Ntot = 0;
            obj.I = zeros(Ntotal, 1);
            obj.J = zeros(Ntotal, 1);
            obj.level = zeros(Ntotal, 1);
            obj.CoarseFactor = zeros(Ntotal, 1);
            obj.CellIndex = zeros(Ntotal, 1);
            obj.Fathers = zeros(Ntotal,1);
            obj.Children = cell(Ntotal, 1);
            obj.GrandChildren = cell(Ntotal, 1);
            %obj.Centers = zeros(Ntotal);
        end
    end
end