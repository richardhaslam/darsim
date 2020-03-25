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
        function BuildFamily(obj, CoarseGrid, FineGrid, CF, currentLevel, maxLevel) 
            % Assign all the Children and GrandChildren
            for L = 1: currentLevel
                for c = 1:CoarseGrid.N
                    if L == 1
                        %% Coordinates of fine cells contained in the coarse block for children
                        Imin = CoarseGrid.I(c,2) - floor(CF(1)/2); Imin = max(Imin,1);
                        Imax = CoarseGrid.I(c,2) + floor(CF(1)/2); Imax = min(Imax,FineGrid.Nx);
                        Jmin = CoarseGrid.J(c,2) - floor(CF(2)/2); Jmin = max(Jmin,1);
                        Jmax = CoarseGrid.J(c,2) + floor(CF(2)/2); Jmax = min(Jmax,FineGrid.Ny);
                        Kmin = CoarseGrid.K(c,2) - floor(CF(3)/2); Kmin = max(Kmin,1);
                        Kmax = CoarseGrid.K(c,2) + floor(CF(3)/2); Kmax = min(Kmax,FineGrid.Nz);
                        % indeces of the fine cells
                        i = Imin:Imax;
                        j = Jmin:Jmax;
                        k = Kmin:Kmax;
                        [p,q, v] = meshgrid(i, j, k);
                        pairs = [p(:), q(:), v(:)];
                        indexes = pairs(:,1) + (pairs(:,2)-1)*FineGrid.Nx + (pairs(:,3) - 1)*FineGrid.Nx*FineGrid.Ny;
                        CoarseGrid.Children{c,L}(1:length(indexes)) = sort(indexes);
                    else
                        CoarseGrid.Children{c,L} = sort([FineGrid.Children{CoarseGrid.Children{c,L-1},currentLevel}]);
                    end
                end
            end
            for L = currentLevel+1 : maxLevel
            end
        end
        function AssignFathersandVerteces(obj, FineGrid, CoarseGrid, maxLevel)
            % Fine Grid
            FineGrid.Fathers = zeros(FineGrid.N, maxLevel);
            FineGrid.Verteces = zeros(FineGrid.N, maxLevel);
            Children = cell(CoarseGrid(1).N,maxLevel);
            for L=1:maxLevel
                Nc = CoarseGrid(L).N;
                for c = 1:Nc
                    Children{c,L} = CoarseGrid(L).Children{c,1};
                    FineGrid.Fathers(Children{c,L}, L) = c;
                    Imin = CoarseGrid(L).I(c, 2) - floor(CoarseGrid(L).CoarseFactor(1)/2); Imin = max(Imin,1);
                    Imax = CoarseGrid(L).I(c, 2) + floor(CoarseGrid(L).CoarseFactor(1)/2); Imax = min(Imax,FineGrid.Nx);
                    Jmin = CoarseGrid(L).J(c, 2) - floor(CoarseGrid(L).CoarseFactor(2)/2); Jmin = max(Jmin,1);
                    Jmax = CoarseGrid(L).J(c, 2) + floor(CoarseGrid(L).CoarseFactor(2)/2); Jmax = min(Jmax,FineGrid.Ny);
                    Kmin = CoarseGrid(L).K(c, 2) - floor(CoarseGrid(L).CoarseFactor(3)/2); Kmin = max(Kmin,1);
                    Kmax = CoarseGrid(L).K(c, 2) + floor(CoarseGrid(L).CoarseFactor(3)/2); Kmax = min(Kmax,FineGrid.Nz);
                    fc = CoarseGrid(L).I(c, 2) + (CoarseGrid(L).J(c, 2) - 1)*FineGrid.Nx + (CoarseGrid(L).K(c, 2) - 1)*FineGrid.Nx*FineGrid.Ny;
                    FineGrid.Verteces(fc, L) = 1;
                    %Scan fine cells inside c
                    for I = Imin:Imax
                        for J = Jmin:Jmax
                            for K = Kmin:Kmax
                                f = I + (J-1)*FineGrid.Nx + (K-1)*FineGrid.Nx*FineGrid.Ny;
                                FineGrid.Fathers(f, L) = c;
                            end
                        end
                    end
                end
            end
            % Coarse Grids
            for x = 1:maxLevel-1
                Nf = CoarseGrid(x).N;
                CoarseGrid(x).Fathers = zeros(Nf, maxLevel);
                CoarseGrid(x).Verteces = zeros(Nf, maxLevel);
                CoarseGrid(x).Fathers(:,x) = ones(Nf, 1) .* [1:CoarseGrid(x).N]' ;
                CoarseGrid(x).Verteces(:,x) = ones(Nf, 1);
                for y = x+1:maxLevel
                    Nc = CoarseGrid(y).N;
                    for c = 1:Nc
                        CoarseGrid(x).Fathers(CoarseGrid(y).Children{c, :}, y) = c;
                        for child = CoarseGrid(y).Children{c, :}
                            if CoarseGrid(y).I(c, 2) == CoarseGrid(x).I(child, 2) && CoarseGrid(y).J(c, 2) == CoarseGrid(x).J(child, 2) && CoarseGrid(y).K(c, 2) == CoarseGrid(x).K(child, 2)
                                CoarseGrid(x).Verteces(child, y) = 1;
                            end
                        end
                    end
                end
            end
            CoarseGrid(maxLevel).Fathers = ones(CoarseGrid(maxLevel).N, maxLevel) .* [zeros(CoarseGrid(maxLevel).N, maxLevel - 1), (1:CoarseGrid(maxLevel).N)'];
            CoarseGrid(maxLevel).Verteces = ones(CoarseGrid(maxLevel).N, maxLevel);
        end
    end
end