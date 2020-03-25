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
                %% Coordinates of fine cells contained in the coarse block for children
                if ~CoarseGrid.Vertex_On_Corner
                    I = CoarseGrid.I(c, 1);
                    J = CoarseGrid.J(c, 1);
                    K = CoarseGrid.K(c, 1);
                    Imin = (I - 1) * CF(1) + 1;
                    Imax = Imin + CF(1) - 1;
                    Jmin = (J - 1) * CF(2) + 1;
                    Jmax = Jmin + CF(2) - 1;
                    Kmin = (K - 1) * CF(3) + 1;
                    Kmax = Kmin + CF(3) - 1;
                else
                    Imin = CoarseGrid.I(c,2) - floor(CF(1)/2); Imin = max(Imin,1);
                    Imax = CoarseGrid.I(c,2) + floor(CF(1)/2); Imax = min(Imax,FineGrid.Nx);
                    Jmin = CoarseGrid.J(c,2) - floor(CF(2)/2); Jmin = max(Jmin,1);
                    Jmax = CoarseGrid.J(c,2) + floor(CF(2)/2); Jmax = min(Jmax,FineGrid.Ny);
                    Kmin = CoarseGrid.K(c,2) - floor(CF(3)/2); Kmin = max(Kmin,1);
                    Kmax = CoarseGrid.K(c,2) + floor(CF(3)/2); Kmax = min(Kmax,FineGrid.Nz);
                end
                
                % indeces of the fine cells
                i = Imin:Imax;
                j = Jmin:Jmax;
                k = Kmin:Kmax;
                [p,q, v] = meshgrid(i, j, k);
                pairs = [p(:), q(:), v(:)];
                indexes = pairs(:,1) + (pairs(:,2)-1)*FineGrid.Nx + (pairs(:,3) - 1)*FineGrid.Nx*FineGrid.Ny;
                CoarseGrid.Children{c,1}(1:length(indexes)) = sort(indexes);
                
                %% Coordinates of fine cells contained in the coarse block for grandchildren
                if level > 1
                    if ~CoarseGrid.Vertex_On_Corner
                        FineFineGrid_Nx = CoarseGrid.Nx * CF(1)^level;
                        FineFineGrid_Ny = CoarseGrid.Ny * CF(2)^level;
                        FineFineGrid_Nz = CoarseGrid.Nz * CF(3)^level;
                        Imin = CoarseGrid.I(c,2) - floor((CF(1)^level - 1)/2);
                        Imax = CoarseGrid.I(c,2) + ceil ((CF(1)^level - 1)/2);
                        Jmin = CoarseGrid.J(c,2) - floor((CF(2)^level - 1)/2);
                        Jmax = CoarseGrid.J(c,2) + ceil ((CF(2)^level - 1)/2);
                        Kmin = CoarseGrid.K(c,2) - floor((CF(3)^level - 1)/2);
                        Kmax = CoarseGrid.K(c,2) + ceil ((CF(3)^level - 1)/2);
                    else
                        FineFineGrid_Nx = (CoarseGrid.Nx-1) * CF(1)^level + 1;
                        FineFineGrid_Ny = (CoarseGrid.Ny-1) * CF(2)^level + 1;
                        FineFineGrid_Nz = (CoarseGrid.Nz-1) * CF(3)^level + 1;
                        Imin = ((CoarseGrid.I(c, 2)-1)/CF(1))+1 - floor(CF(1)^2/2); Imin = max(Imin,1);
                        Imax = ((CoarseGrid.I(c, 2)-1)/CF(1))+1 + floor(CF(1)^2/2); Imax = min(Imax,FineFineGrid_Nx);
                        Jmin = ((CoarseGrid.J(c, 2)-1)/CF(2))+1 - floor(CF(2)^2/2); Jmin = max(Jmin,1);
                        Jmax = ((CoarseGrid.J(c, 2)-1)/CF(2))+1 + floor(CF(2)^2/2); Jmax = min(Jmax,FineFineGrid_Ny);
                        Kmin = ((CoarseGrid.K(c, 2)-1)/CF(3))+1 - floor(CF(3)^2/2); Kmin = max(Kmin,1);
                        Kmax = ((CoarseGrid.K(c, 2)-1)/CF(3))+1 + floor(CF(3)^2/2); Kmax = min(Kmax,FineFineGrid_Nz);
                    end
                    
                    % indeces of the fine cells
                    i = Imin:Imax;
                    j = Jmin:Jmax;
                    k = Kmin:Kmax;
                    [p, q, v] = meshgrid(i, j, k);
                    pairs = [p(:), q(:), v(:)];
                    indexes = pairs(:,1) + (pairs(:,2)-1) * FineFineGrid_Nx + (pairs(:,3)-1) * FineFineGrid_Nx * FineFineGrid_Ny;
                    CoarseGrid.GrandChildren{c,1}(1:length(indexes)) = sort(indexes);
                else % if level == 1
                    CoarseGrid.GrandChildren{c,1} = CoarseGrid.Children{c,1};
                end
            end
        end
        function AssignFathersandVerteces(obj, FineGrid, CoarseGrid, maxLevel)
            % Fine Grid
            FineGrid.Fathers = zeros(FineGrid.N, maxLevel);
            FineGrid.Verteces = zeros(FineGrid.N, maxLevel);
            for i=1:maxLevel
                Nc = CoarseGrid(i).N;
                for c = 1:Nc
                    Imin = CoarseGrid(i).I(c, 2) - floor(CoarseGrid(i).CoarseFactor(1)/2); Imin = max(Imin,1);
                    Imax = CoarseGrid(i).I(c, 2) + floor(CoarseGrid(i).CoarseFactor(1)/2); Imax = min(Imax,FineGrid.Nx);
                    Jmin = CoarseGrid(i).J(c, 2) - floor(CoarseGrid(i).CoarseFactor(2)/2); Jmin = max(Jmin,1);
                    Jmax = CoarseGrid(i).J(c, 2) + floor(CoarseGrid(i).CoarseFactor(2)/2); Jmax = min(Jmax,FineGrid.Ny);
                    Kmin = CoarseGrid(i).K(c, 2) - floor(CoarseGrid(i).CoarseFactor(3)/2); Kmin = max(Kmin,1);
                    Kmax = CoarseGrid(i).K(c, 2) + floor(CoarseGrid(i).CoarseFactor(3)/2); Kmax = min(Kmax,FineGrid.Nz);
                    fc = CoarseGrid(i).I(c, 2) + (CoarseGrid(i).J(c, 2) - 1)*FineGrid.Nx + (CoarseGrid(i).K(c, 2) - 1)*FineGrid.Nx*FineGrid.Ny;
                    FineGrid.Verteces(fc, i) = 1;
                    %Scan fine cells inside c
                    for I = Imin:Imax
                        for J = Jmin:Jmax
                            for K = Kmin:Kmax
                                f = I + (J-1)*FineGrid.Nx + (K-1)*FineGrid.Nx*FineGrid.Ny;
                                FineGrid.Fathers(f, i) = c;
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
    