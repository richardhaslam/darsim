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
            %obj.CoarseFactor(obj.Ntot, 1) = [obj.CoarseFactor(1:Nf, 1); zeros(Nx, 1)];
            %obj.CoarseFactor(obj.Ntot, 2) = [obj.CoarseFactor(1:Nf, 2); zeros(Nx, 1)];
            obj.CellIndex = [obj.CellIndex(1:Nf); zeros(Nx, 1)];
            obj.Fathers = [obj.Fathers(1:Nf,:); zeros(Nx, 2)];
            obj.Verteces = [obj.Verteces(1:Nf,:); zeros(Nx, 2)];
            obj.level = [obj.level(1:Nf); ones(Nx, 1)*obj.MaxLevel];
            
            % Add the new cells
            h = Nf + 1;
            for CoarseNode = Nf+1 : Nf+Nc;
                FineNodes = obj.Children{CoarseNode};
                obj.I(h:h + 8) = FineGrid.I(FineNodes);
                obj.J(h:h + 8) = FineGrid.J(FineNodes);
                %obj.CoarseFactor(CoarseNode:CoarseNode + 9, 1) = FineGrid.CoarseFactor(1);
                %obj.CoarseFactor(CoarseNode:CoarseNode + 9, 2) = FineGrid.CoarseFactor(2);
                obj.CellIndex(h:h + 8) = FineNodes;
                obj.Fathers(h:h + 8,:) = FineGrid.Fathers(FineNodes,:);
                obj.Verteces(h:h + 8,:) = FineGrid.Verteces(FineNodes,:);
                h = h + 9;
            end
        end
    end
end