%  ADM discretization model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 12 July 2016
%Last modified: 26 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef ADM_Discretization_model < Discretization_model
    properties
        maxLevel
        Coarsening
        CoarseGrid
        ADMGrid
        ADMGridSelector
        OperatorsHandler
    end
    methods
        function obj = ADM_Discretization_model(nx, ny, nz, maxlevel, coarsening)
            obj@Discretization_model(nx, ny, nz);
            obj.CoarseGrid = coarse_grid();
            obj.Coarsening = coarsening;
            obj.maxLevel = maxlevel;
        end
        function AddOperatorsHandler(obj, operatorshandler)
            obj.OperatorsHandler = operatorshandler;
        end
        function AddADMGridSelector(obj, gridselector)
            obj.ADMGridSelector = gridselector;
        end
        function Initialize(obj, ProductionSystem)
            disp('Algebraic Dynamic Multilevel (ADM) method run')
            obj.ReservoirGrid.Initialize(ProductionSystem.Reservoir);
            % Perforated cells
            obj.DefinePerforatedCells(ProductionSystem.Wells);
            
            % Construct Coarse Grids
            obj.ConstructCoarseGrids(ProductionSystem.Wells.Inj, ProductionSystem.Wells.Prod);
            
            %% Pressure interpolators
            disp('Static operators - start computation');
            obj.OperatorsHandler.BuildStaticOperators(obj.CoarseGrid, obj.ReservoirGrid, obj.maxLevel, ProductionSystem.Reservoir.K);
            disp('Static operators - end')
            disp(char(2));
        end
        function ConstructCoarseGrids(obj, Inj, Prod)
            % Add coordinates to fine-scale grid (useful for ADM)
            obj.ReservoirGrid.AddCoordinates();
            
            % Construct all coarse grids
            for i=1:obj.maxLevel
                obj.CoarseGrid(i).CoarseFactor = obj.Coarsening(i,:);
                obj.CoarseGrid(i).BuildCoarseGrid(obj.ReservoirGrid);
            end
            
            % Assign fathers and verteces 
            obj.AssignFathers();
            
            % Flag coarse blocks with wells
            obj.CoarseWells(Inj, Prod);
            obj.ADMGridSelector.NoWellsCoarseCells = ones(obj.CoarseGrid(1).N, 1);
            Nc1 = obj.CoarseGrid(1).N;
            if obj.maxLevel > 1
               for i = 1:Nc1
                   if obj.CoarseGrid(2).Wells(obj.CoarseGrid(1).Father(i,2)) == 1
                       obj.ADMGridSelector.NoWellsCoarseCells(i) = 0;
                   end
               end
            else
               for i =1:Nc1
                   if obj.CoarseGrid(1).Wells(i) == 1
                       obj.ADMGridSelector.NoWellsCoarseCells(i) = 0;
                   end
               end
            end
        end
        function AssignFathers(obj)
            % Fine Grid
            obj.ReservoirGrid.Father = zeros(obj.ReservoirGrid.N, obj.maxLevel);
            obj.ReservoirGrid.Centers = zeros(obj.ReservoirGrid.N, obj.maxLevel);
            for i=1:obj.maxLevel
                Nc = obj.CoarseGrid(i).Nx*obj.CoarseGrid(i).Ny;
                for c = 1:Nc
                    Imin = obj.CoarseGrid(i).I(c) - floor((obj.CoarseGrid(i).CoarseFactor(1) - 1)/2);
                    Imax = obj.CoarseGrid(i).I(c) + ceil((obj.CoarseGrid(i).CoarseFactor(1) - 1)/2);
                    Jmin = obj.CoarseGrid(i).J(c) - floor((obj.CoarseGrid(i).CoarseFactor(2) - 1)/2);
                    Jmax = obj.CoarseGrid(i).J(c) + ceil((obj.CoarseGrid(i).CoarseFactor(2) - 1)/2);
                    fc = obj.CoarseGrid(i).I(c) + (obj.CoarseGrid(i).J(c) - 1)*obj.ReservoirGrid.Nx;
                    obj.ReservoirGrid.Centers(fc,i) = 1;
                    %Scan fine cells inside c
                    for I = Imin:Imax
                        for J = Jmin:Jmax
                            f = I + (J-1)*obj.ReservoirGrid.Nx;
                            obj.ReservoirGrid.Father(f,i) = c;
                        end
                    end
                end
                obj.CoarseGrid(i).Centers = ones(obj.CoarseGrid(i).Nx*obj.CoarseGrid(i).Ny, 1);
            end
            % Coarse Grids
            for x = 1:obj.maxLevel-1
                Nf = obj.CoarseGrid(x).Nx*obj.CoarseGrid(x).Ny;
                obj.CoarseGrid(x).Father = zeros(Nf, obj.maxLevel);
                obj.CoarseGrid(x).Centers = zeros(Nf, obj.maxLevel);
                %obj.CoarseGrid(x).Centers(:, x) = ones(Nf, 1);
                for y = x+1:obj.maxLevel
                    Nc = obj.CoarseGrid(y).Nx*obj.CoarseGrid(y).Ny;
                    for c = 1:Nc
                        Imin = obj.CoarseGrid(y).I(c) - floor((obj.CoarseGrid(y).CoarseFactor(1) - 1)/2);
                        Imax = obj.CoarseGrid(y).I(c) + ceil((obj.CoarseGrid(y).CoarseFactor(1) - 1)/2);
                        Jmin = obj.CoarseGrid(y).J(c) - floor((obj.CoarseGrid(y).CoarseFactor(2) - 1)/2);
                        Jmax = obj.CoarseGrid(y).J(c) + ceil((obj.CoarseGrid(y).CoarseFactor(2) - 1)/2);
                        for l = 1:Nf
                            if ((obj.CoarseGrid(x).I(l) <= Imax && obj.CoarseGrid(x).I(l) >= Imin) && ...
                                    (obj.CoarseGrid(x).J(l) <= Jmax && obj.CoarseGrid(x).J(l) >= Jmin))
                                obj.CoarseGrid(x).Father(l,y) = c;
                            end
                            if (obj.CoarseGrid(x).I(l) == obj.CoarseGrid(y).I(c) && obj.CoarseGrid(x).J(l) == obj.CoarseGrid(y).J(c))
                                obj.CoarseGrid(x).Centers(l, y) = 1;
                            end
                        end
                    end
                end
            end
            obj.CoarseGrid(obj.maxLevel).Father = zeros(obj.CoarseGrid(obj.maxLevel).N, obj.maxLevel);
            obj.CoarseGrid(obj.maxLevel).Centers = zeros(obj.CoarseGrid(obj.maxLevel).N, obj.maxLevel);
        end
        function CoarseWells(obj, Inj, Prod)
            for i=1:length(Inj)
                % Flag coarse Nodes with wells
                I = Inj(i).Cells;
                for x = 1:obj.maxLevel
                    obj.CoarseGrid(x).Wells(obj.ReservoirGrid.Father(I,x)) = 1;
                end
            end
            for i =1:length(Prod)
                P = Prod(i).Cells;
                for x = 1:obj.maxLevel
                    obj.CoarseGrid(x).Wells(obj.ReservoirGrid.Father(P,x)) = 1;
                end
            end
        end      
    end
end