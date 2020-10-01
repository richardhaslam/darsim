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
        function BuildFamily(obj, CoarseGrid, FineGrid, CF, currentLevel) 
            % Assign all the Children and GrandChildren
            for L = 1 : currentLevel
                for c = 1:CoarseGrid.N
                    if L == 1
                        %% Coordinates of fine cells contained in the coarse block for children
                        Imin = CoarseGrid.I(c,2) - floor(CF(1)/2); Imin = max(Imin,1);
                        Imax = CoarseGrid.I(c,2) + floor(CF(1)/2); Imax = min(Imax,FineGrid(end).Nx);
                        Jmin = CoarseGrid.J(c,2) - floor(CF(2)/2); Jmin = max(Jmin,1);
                        Jmax = CoarseGrid.J(c,2) + floor(CF(2)/2); Jmax = min(Jmax,FineGrid(end).Ny);
                        Kmin = CoarseGrid.K(c,2) - floor(CF(3)/2); Kmin = max(Kmin,1);
                        Kmax = CoarseGrid.K(c,2) + floor(CF(3)/2); Kmax = min(Kmax,FineGrid(end).Nz);
                        % indeces of the fine cells
                        i = Imin:Imax;
                        j = Jmin:Jmax;
                        k = Kmin:Kmax;
                        [p,q, v] = meshgrid(i, j, k);
                        pairs = [p(:), q(:), v(:)];
                        indexes = pairs(:,1) + (pairs(:,2)-1)*FineGrid(end).Nx + (pairs(:,3) - 1)*FineGrid(end).Nx*FineGrid(end).Ny;
                        CoarseGrid.Children{c,L}(1:length(indexes)) = sort(indexes);
                    else
                        CoarseGrid.Children{c,L} = sort([FineGrid(end-L+2).Children{CoarseGrid.Children{c,L-1},1}]);
                    end
                end
            end
        end
        function AssignFathersandVerteces(obj, Grid, maxLevel)
            for L=1:maxLevel
                Grid(L).Fathers = zeros(Grid(L).N, maxLevel-L+1);
                Grid(L).Verteces = zeros(Grid(L).N, maxLevel-L+1);
                for LL = 1 : maxLevel-L+1
                    Nc = Grid(LL+L).N;
                    for c = 1:Nc
                        % Assigning Fathers
                        Grid(L).Fathers(Grid(LL+L).Children{c,LL}, LL) = c;
                        % Assigning Verteces
                        CF = zeros(3,1);
                        if Grid(LL+L).Vertex_On_Corner
                            CF(1) = (Grid(L).Nx-1)/(Grid(LL+L-1).Nx-1);
                            CF(2) = (Grid(L).Ny-1)/(Grid(LL+L-1).Ny-1);
                            CF(3) = (Grid(L).Nz-1)/(Grid(LL+L-1).Nz-1);
                            CF(isnan(CF))=1;
                        else
                            CF(1) = Grid(L).Nx/Grid(LL+L-1).Nx;
                            CF(2) = Grid(L).Ny/Grid(LL+L-1).Ny;
                            CF(3) = Grid(L).Nz/Grid(LL+L-1).Nz;
                        end
                        Imin = Grid(LL+L).I(c,LL+1) - floor(CF(1)/2); Imin = max(Imin,1);
                        Imax = Grid(LL+L).I(c,LL+1) + floor(CF(1)/2); Imax = min(Imax,Grid(L).Nx);
                        Jmin = Grid(LL+L).J(c,LL+1) - floor(CF(2)/2); Jmin = max(Jmin,1);
                        Jmax = Grid(LL+L).J(c,LL+1) + floor(CF(2)/2); Jmax = min(Jmax,Grid(L).Ny);
                        Kmin = Grid(LL+L).K(c,LL+1) - floor(CF(3)/2); Kmin = max(Kmin,1);
                        Kmax = Grid(LL+L).K(c,LL+1) + floor(CF(3)/2); Kmax = min(Kmax,Grid(L).Nz);
                        i = Imin:Imax;
                        j = Jmin:Jmax;
                        k = Kmin:Kmax;
                        [p,q, v] = meshgrid(i, j, k);
                        pairs = [p(:), q(:), v(:)];
                        VertexIndex = sort( pairs(:,1) + (pairs(:,2)-1)*Grid(L).Nx + (pairs(:,3) - 1)*Grid(L).Nx*Grid(L).Ny );
                        Grid(L).Verteces(VertexIndex,LL) = 1; 
                    end
                end
            end
        end
    end
end