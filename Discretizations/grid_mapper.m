%  Grid mapper 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 5 September 2016
%Last modified: 5 September 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef grid_mapper < handle
    properties
    end
    methods
        function BuildFamily(obj, CoarseGrid, FineGrid, CF, level) 
        % Assign Children and GrandChildren 
          for c = 1:CoarseGrid.N
              %coordinates of fine cells contained in the coarse block
              J = ceil(c/CoarseGrid.Ny);
              I = c - CoarseGrid.Nx*(J - 1);
              Imin = (I - 1) * CF(2) + 1;
              Imax = Imin + CF(1)-1;
              Jmin = (J - 1) * CF(1) + 1;
              Jmax = Jmin + CF(2)-1;
              %indexes of the fine cells
              i = Imin:Imax;
              j = Jmin:Jmax;
              [p,q] = meshgrid(i, j);
              pairs = [p(:), q(:)];
              indexes = pairs(:,1) + (pairs(:,2)-1)*FineGrid.Nx;
              CoarseGrid.Children(c, :) = indexes;
              
              %coordinates of fine cells contained in the coarse block
              Imin = CoarseGrid.I(c) - floor((CF(1)^level - 1)/2);
              Imax = CoarseGrid.I(c) + ceil((CF(1)^level - 1)/2);
              Jmin = CoarseGrid.J(c) - floor((CF(2)^level - 1)/2);
              Jmax = CoarseGrid.J(c) + ceil((CF(2)^level - 1)/2);
              i = Imin:Imax;
              j = Jmin:Jmax;
              [p,q] = meshgrid(i, j);
              pairs = [p(:), q(:)];
              %indexes of the fine cells
              indexes = pairs(:,1) + (pairs(:,2)-1)*CoarseGrid.Nx*CF(1)^level;
              CoarseGrid.GrandChildren(c, :) = indexes;
          end          
        end
        function AssignFathersandVerteces(obj, FineGrid, CoarseGrid, maxLevel)
            % Fine Grid
            FineGrid.Fathers = zeros(FineGrid.N, maxLevel);
            FineGrid.Verteces = zeros(FineGrid.N, maxLevel);
            for i=1:maxLevel
                Nc = CoarseGrid(i).Nx*CoarseGrid(i).Ny;
                for c = 1:Nc
                    Imin = CoarseGrid(i).I(c) - floor((CoarseGrid(i).CoarseFactor(1) - 1)/2);
                    Imax = CoarseGrid(i).I(c) + ceil((CoarseGrid(i).CoarseFactor(1) - 1)/2);
                    Jmin = CoarseGrid(i).J(c) - floor((CoarseGrid(i).CoarseFactor(2) - 1)/2);
                    Jmax = CoarseGrid(i).J(c) + ceil((CoarseGrid(i).CoarseFactor(2) - 1)/2);
                    fc = CoarseGrid(i).I(c) + (CoarseGrid(i).J(c) - 1)*FineGrid.Nx;
                    FineGrid.Verteces(fc, i) = 1;
                    %Scan fine cells inside c
                    for I = Imin:Imax
                        for J = Jmin:Jmax
                            f = I + (J-1)*FineGrid.Nx;
                            FineGrid.Fathers(f,i) = c;
                        end
                    end
                end
                CoarseGrid(i).Verteces = ones(CoarseGrid(i).Nx*CoarseGrid(i).Ny, 1);
            end
            % Coarse Grids
            for x = 1:maxLevel-1
                Nf = CoarseGrid(x).N;
                CoarseGrid(x).Fathers = zeros(Nf, maxLevel);
                CoarseGrid(x).Verteces = zeros(Nf, maxLevel);
                for y = x+1:maxLevel
                    Nc = CoarseGrid(y).N;
                    for c = 1:Nc
                        Imin = CoarseGrid(y).I(c) - floor((CoarseGrid(y).CoarseFactor(1) - 1)/2);
                        Imax = CoarseGrid(y).I(c) + ceil((CoarseGrid(y).CoarseFactor(1) - 1)/2);
                        Jmin = CoarseGrid(y).J(c) - floor((CoarseGrid(y).CoarseFactor(2) - 1)/2);
                        Jmax = CoarseGrid(y).J(c) + ceil((CoarseGrid(y).CoarseFactor(2) - 1)/2);
                        for l = 1:Nf
                            if ((CoarseGrid(x).I(l) <= Imax && CoarseGrid(x).I(l) >= Imin) && ...
                                    (CoarseGrid(x).J(l) <= Jmax && CoarseGrid(x).J(l) >= Jmin))
                                CoarseGrid(x).Fathers(l,y) = c;
                            end
                            if (CoarseGrid(x).I(l) == CoarseGrid(y).I(c) && CoarseGrid(x).J(l) == CoarseGrid(y).J(c))
                                CoarseGrid(x).Verteces(l, y) = 1;
                            end
                        end
                    end
                end
            end
            CoarseGrid(maxLevel).Fathers = zeros(CoarseGrid(maxLevel).N, maxLevel);
            CoarseGrid(maxLevel).Verteces = zeros(CoarseGrid(maxLevel).N, maxLevel);
        end
    end
end
    