%  dynamic grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 12 July 2016
%Last modified: 12 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef adm_grid < grid_darsim
    properties
        Ntot
        I
        J
        K
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
            obj.K = zeros(Ntotal, 1);
            obj.level = zeros(Ntotal, 1);
            obj.CoarseFactor = zeros(Ntotal, 1);
            obj.CellIndex = zeros(Ntotal, 1);
            obj.Fathers = zeros(Ntotal, maxlevel);
            obj.Children = cell(Ntotal, 1);
            obj.GrandChildren = cell(Ntotal, 1);
            obj.Verteces = zeros(Ntotal, maxlevel);
        end
        function Update(obj, Nx, Nf, Nc, FineGrid)
            obj.N(obj.MaxLevel + 1) = 0;
            obj.MaxLevel = obj.MaxLevel - 1;
            obj.N(obj.MaxLevel + 1) = obj.N(obj.MaxLevel + 1) + Nx;
            obj.Ntot = Nf + Nx;
            
            obj.I = [obj.I(1:Nf); zeros(Nx, 1)];
            obj.J = [obj.J(1:Nf); zeros(Nx, 1)];
            obj.CellIndex = [obj.CellIndex(1:Nf); zeros(Nx, 1)];
            obj.Fathers = [obj.Fathers(1:Nf,:); zeros(Nx, 2)];
            obj.Verteces = [obj.Verteces(1:Nf,:); zeros(Nx, 2)];
            obj.level = [obj.level(1:Nf); ones(Nx, 1)*obj.MaxLevel];
            
            % Add the new cells
            h = Nf + 1;
            for CoarseNode = Nf+1 : Nf+Nc;
                FineNodes = obj.Children{CoarseNode};
                n_children = length(FineNodes);
                obj.I(h:h + n_children-1) = FineGrid.I(FineNodes);
                obj.J(h:h + n_children-1) = FineGrid.J(FineNodes);
                obj.CellIndex(h:h + n_children-1) = FineNodes;
                obj.Fathers(h:h + n_children-1,:) = FineGrid.Fathers(FineNodes,:);
                obj.Verteces(h:h + n_children-1,:) = FineGrid.Verteces(FineNodes,:);
                h = h + n_children;
            end
        end
    end
end