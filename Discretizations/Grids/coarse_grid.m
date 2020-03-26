%  Coarse grid class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini & Mousa HosseiniMehr
%TU Delft
%Created: 26 July 2016
%Last modified: 27 March 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef coarse_grid < grid_darsim
    properties
        Nx
        Ny
        Nz
        I
        J
        K
        Neighbours
        Wells
        DeltaS
        Vertex_On_Corner
    end
    methods
        function obj = coarse_grid() 
        end
        function BuildCoarseGrid(obj, FineGrids, FullCF)            
            %% Construct a coarse Grid given all the coarsening ratios and FineGrids
            FullCF(1,:,2:end+1) = FullCF(1,:,1:end);
            FullCF(1,:,1) = [1,1,1];
            CF = FullCF(:,:,end)./FullCF(:,:,end-1);
            if ~obj.Vertex_On_Corner
                obj.Nx = FineGrids(end).Nx/CF(1);
                obj.Ny = FineGrids(end).Ny/CF(2);
                obj.Nz = FineGrids(end).Nz/CF(3);
            else
                obj.Nx = max( (FineGrids(end).Nx-1)/CF(1)+1 , 1);
                obj.Ny = max( (FineGrids(end).Ny-1)/CF(2)+1 , 1);
                obj.Nz = max( (FineGrids(end).Nz-1)/CF(3)+1 , 1);
            end
            obj.N = obj.Nx*obj.Ny*obj.Nz;
            obj.DeltaS = zeros(obj.N, 1);
            
            % Coordinates of the centres
            obj.I = ones(obj.N, 1+length(FineGrids));
            obj.J = ones(obj.N, 1+length(FineGrids));
            obj.K = ones(obj.N, 1+length(FineGrids));
            
            for L = 1:length(FineGrids)
                FineGrid = FineGrids(end-L+1);
                CF = FullCF(:,:,end)./FullCF(:,:,end-L);
                if ~obj.Vertex_On_Corner
                    JindecesF = ceil(CF(2)/2):CF(2):FineGrid.Ny;
                else
                    JindecesF = 1:CF(2):FineGrid.Ny;
                end
                JindecesC = 1:1:obj.Ny;
                for k = 1:obj.Nz
                    for i=1:obj.Ny
                        a = obj.Nx*(i-1) + (k-1)*obj.Nx*obj.Ny + 1;
                        if ~obj.Vertex_On_Corner
                            obj.I(a:a+obj.Nx-1, 1+L) = ceil(CF(1)/2):CF(1):FineGrid.Nx;
                        else
                            obj.I(a:a+obj.Nx-1, 1+L) = 1:CF(1):FineGrid.Nx;
                        end
                        obj.J(a:a+obj.Nx-1, 1+L) = JindecesF(i)*ones(obj.Nx,1);
                        obj.I(a:a+obj.Nx-1, 1) = 1:1:obj.Nx;
                        obj.J(a:a+obj.Nx-1, 1) = JindecesC(i)*ones(obj.Nx, 1);
                    end
                end
                if ~obj.Vertex_On_Corner
                    KindecesC = ceil(CF(3)/2):CF(3):FineGrid.Nz;
                else
                    KindecesC = 1:CF(3):FineGrid.Nz;
                end
                for i=1:obj.Nz
                    obj.K((i-1)*obj.Nx*obj.Ny+1:i*obj.Nx*obj.Ny, 1+L) = KindecesC(i) * ones(obj.Nx*obj.Ny, 1);
                    obj.K((i-1)*obj.Nx*obj.Ny+1:i*obj.Nx*obj.Ny, 1) = i*ones(obj.Nx*obj.Ny, 1);
                end
            end

            obj.Active = zeros(obj.N, 1);
            obj.Wells = cell(obj.N, 1);
            obj.Fathers = zeros(obj.N, 1);
            
            if obj.Nz == 1 && obj.Ny == 1
                obj.AssignNeighbours1D();
            elseif obj.Nz == 1 && obj.Ny > 1
                obj.AssignNeighbours2D();
            else
                obj.AssignNeighbours();
            end
            
            CF = FullCF(:,:,end)./FullCF(:,:,end-1);
            obj.AddGridCoordinates(FineGrids(end),CF); 
        end
        function AddGridCoordinates(obj, FineGrid, CF)
            % Computes coordinates of corners of the coarse grid
            if ~isempty(FineGrid.GridCoords)
                if ~obj.Vertex_On_Corner
                    kList = 1 : obj.Nz+1;
                    jList = 1 : obj.Ny+1;
                else
                    kList = [0, ceil(CF(3)/2) + CF(3)*(0:obj.Nz-2), FineGrid.Nz] + 1;
                    jList = [0, ceil(CF(2)/2) + CF(2)*(0:obj.Ny-2), FineGrid.Ny] + 1;
                end

                ThreeDim = (size(FineGrid.GridCoords,1) == (FineGrid.Nx+1)*(FineGrid.Ny+1)*(FineGrid.Nz+1));
                if ~ThreeDim
                    kLoopNum = 1;
                else
                    kLoopNum = length(kList);
                end

                for m = 1 : kLoopNum
                    k = kList(m);
                    for n = 1 : length(jList)
                        j = jList(n);
                        if ~obj.Vertex_On_Corner
                            fs_index = CF(3)*(FineGrid.Ny+1)*(FineGrid.Nx+1)*(k-1) + ...
                                       CF(2)*(FineGrid.Nx+1)*(j-1) + ...
                                      (CF(1)*(0:obj.Nx) ) + 1;
                        else
                            fs_index = (k-1)*(FineGrid.Nx+1)*(FineGrid.Ny+1) + (j-1)*(FineGrid.Nx+1) + ...
                                       [0, ceil(CF(1)/2) + CF(1)*(0:obj.Nx-2), FineGrid.Nx] + 1;
                        end
                        
                        Start = (m-1)*(obj.Nx+1)*(obj.Ny+1) + (n-1)*(obj.Nx+1) + 1 ;
                        End   = (m-1)*(obj.Nx+1)*(obj.Ny+1) + (n  )*(obj.Nx+1)     ;
                        obj.GridCoords(Start:End,:) = FineGrid.GridCoords(fs_index,:);
                    end
                end
            end
        end
        function AssignNeighbours(obj)
            % Let s do the 8 corners separetely
            % 1
            i = 1;
            j = 1;
            for k=1:obj.Nz
                g = i + (j-1)*obj.Nx + (k-1)*obj.Nx*obj.Ny;
                if k == 1
                    obj.Neighbours(g).indexes = [g+1, g+obj.Nx, g+obj.Nx*obj.Ny];
                elseif k == obj.Nz
                    obj.Neighbours(g).indexes = [g+1, g+obj.Nx, g-obj.Nx*obj.Ny];
                else
                    obj.Neighbours(g).indexes = [g+1, g+obj.Nx, g-obj.Nx*obj.Ny, g+obj.Nx*obj.Ny];
                end
            end
            % 2.
            i = obj.Nx;
            j = 1;
            for k=1:obj.Nz
                g = i + (j-1)*obj.Nx + (k-1)*obj.Nx*obj.Ny;
                if k == 1
                    obj.Neighbours(g).indexes = [g-1, g+obj.Nx, g+obj.Nx*obj.Ny];
                elseif k == obj.Nz
                    obj.Neighbours(g).indexes = [g-1, g+obj.Nx, g-obj.Nx*obj.Ny];
                else
                    obj.Neighbours(g).indexes = [g-1, g+obj.Nx, g-obj.Nx*obj.Ny, g+obj.Nx*obj.Ny];
                end
            end
            % 3.
            i = 1;
            j = obj.Ny;
            for k=1:obj.Nz
                g = i + (j-1)*obj.Nx + (k-1)*obj.Nx*obj.Ny;
                if k == 1
                    obj.Neighbours(g).indexes = [g+1, g-obj.Nx, g+obj.Nx*obj.Ny];
                elseif k == obj.Nz
                    obj.Neighbours(g).indexes = [g+1, g-obj.Nx, g-obj.Nx*obj.Ny];
                else
                    obj.Neighbours(g).indexes = [g+1, g-obj.Nx, g-obj.Nx*obj.Ny, g+obj.Nx*obj.Ny];
                end
            end
            % 4.
            i = obj.Nx;
            j = obj.Ny;
            for k=1:obj.Nz
                g = i + (j-1)*obj.Nx + (k-1)*obj.Nx*obj.Ny;
                if k == 1
                    obj.Neighbours(g).indexes = [g-1, g-obj.Nx, g+obj.Nx*obj.Ny];
                elseif k == obj.Nz
                    obj.Neighbours(g).indexes = [g-1, g-obj.Nx, g-obj.Nx*obj.Ny];
                else
                    obj.Neighbours(g).indexes = [g-1, g-obj.Nx, g-obj.Nx*obj.Ny, g+obj.Nx*obj.Ny];
                end
            end
            % 5
            j = 1;
            k = 1;
            for i=1:obj.Nx
                g = i + (j-1)*obj.Nx + (k-1)*obj.Nx*obj.Ny;
                if i==1
                    obj.Neighbours(g).indexes = [g+1, g+obj.Nx, g+obj.Nx*obj.Ny];
                elseif i==obj.Nx
                    obj.Neighbours(g).indexes = [g-1, g+obj.Nx, g+obj.Nx*obj.Ny];
                else
                    obj.Neighbours(g).indexes = [g-1, g+1, g+obj.Nx, g+obj.Nx*obj.Ny];
                end
            end
            % 6
            j = obj.Ny;
            k = 1;
            for i=1:obj.Nx
                g = i + (j-1)*obj.Nx + (k-1)*obj.Nx*obj.Ny;
                if i==1
                    obj.Neighbours(g).indexes = [g+1, g-obj.Nx, g+obj.Nx*obj.Ny];
                elseif i==obj.Nx
                    obj.Neighbours(g).indexes = [g-1, g-obj.Nx, g+obj.Nx*obj.Ny];
                else
                    obj.Neighbours(g).indexes = [g-1, g+1, g-obj.Nx, g+obj.Nx*obj.Ny];
                end
            end
            % 7
            j = 1;
            k = obj.Nz;
            for i=1:obj.Nx
                g = i + (j-1)*obj.Nx + (k-1)*obj.Nx*obj.Ny;
                if i==1
                    obj.Neighbours(g).indexes = [g+1, g+obj.Nx, g-obj.Nx*obj.Ny];
                elseif i==obj.Nx
                    obj.Neighbours(g).indexes = [g-1, g+obj.Nx, g-obj.Nx*obj.Ny];
                else
                    obj.Neighbours(g).indexes = [g-1, g+1, g+obj.Nx, g-obj.Nx*obj.Ny];
                end
            end
            % 8
            j = obj.Ny;
            k = obj.Nz;
            for i=1:obj.Nx
                g = i + (j-1)*obj.Nx + (k-1)*obj.Nx*obj.Ny;
                if i==1
                    obj.Neighbours(g).indexes = [g+1, g-obj.Nx, g-obj.Nx*obj.Ny];
                elseif i==obj.Nx
                    obj.Neighbours(g).indexes = [g-1, g-obj.Nx, g-obj.Nx*obj.Ny];
                else
                    obj.Neighbours(g).indexes = [g-1, g+1, g-obj.Nx, g-obj.Nx*obj.Ny];
                end
            end
            
            % 9
            i = 1;
            k = 1;
            for j=1:obj.Ny
                g = i + (j-1)*obj.Nx + (k-1)*obj.Nx*obj.Ny;
                if j==1
                    obj.Neighbours(g).indexes = [g+1, g+obj.Nx, g+obj.Nx*obj.Ny];
                elseif j==obj.Ny
                    obj.Neighbours(g).indexes = [g-1, g-obj.Nx, g+obj.Nx*obj.Ny];
                else
                    obj.Neighbours(g).indexes = [g+1, g-obj.Nx, g+obj.Nx, g+obj.Nx*obj.Ny];
                end
            end
            % 10
            i = obj.Nx;
            k = 1;
            for j=1:obj.Ny
                g = i + (j-1)*obj.Nx + (k-1)*obj.Nx*obj.Ny;
                if j==1
                    obj.Neighbours(g).indexes = [g-1, g+obj.Nx, g+obj.Nx*obj.Ny];
                elseif j==obj.Ny
                    obj.Neighbours(g).indexes = [g-1, g-obj.Nx, g+obj.Nx*obj.Ny];
                else
                    obj.Neighbours(g).indexes = [g-1, g-obj.Nx, g+obj.Nx, g+obj.Nx*obj.Ny];
                end
            end
            % 11
            i = 1;
            k = obj.Nz;
            for j=1:obj.Ny
                g = i + (j-1)*obj.Nx + (k-1)*obj.Nx*obj.Ny;
                if j==1
                    obj.Neighbours(g).indexes = [g+1, g+obj.Nx, g-obj.Nx*obj.Ny];
                elseif j==obj.Ny
                    obj.Neighbours(g).indexes = [g+1, g-obj.Nx, g-obj.Nx*obj.Ny];
                else
                    obj.Neighbours(g).indexes = [g+1, g-obj.Nx, g+obj.Nx, g-obj.Nx*obj.Ny];
                end
            end
            % 12
            i = obj.Nx;
            k = obj.Nz;
            for j=1:obj.Ny
                g = i + (j-1)*obj.Nx + (k-1)*obj.Nx*obj.Ny;
                if j==1
                    obj.Neighbours(g).indexes = [g-1, g+obj.Nx, g-obj.Nx*obj.Ny];
                elseif j==obj.Ny
                    obj.Neighbours(g).indexes = [g-1, g-obj.Nx, g-obj.Nx*obj.Ny];
                else
                    obj.Neighbours(g).indexes = [g-1, g-obj.Nx, g+obj.Nx, g-obj.Nx*obj.Ny];
                end
            end

            % Let's do the faces first
            i = 1;
            for k=2:obj.Nz-1
                for j=2:obj.Ny-1
                    g = i + (j-1)*obj.Nx + (k-1)*obj.Nx*obj.Ny;
                    obj.Neighbours(g).indexes = [g+1, g-obj.Nx, g+obj.Nx, g-obj.Nx*obj.Ny, g+obj.Nx*obj.Ny];
                end
            end
            i = obj.Nx;
            for k=2:obj.Nz-1
                for j=2:obj.Ny-1
                    g = i + (j-1)*obj.Nx + (k-1)*obj.Nx*obj.Ny;
                    obj.Neighbours(g).indexes = [g-1, g-obj.Nx, g+obj.Nx, g-obj.Nx*obj.Ny, g + obj.Nx*obj.Ny];
                end
            end
            j = 1;
            for k=2:obj.Nz-1
                for i=2:obj.Nx-1
                    g = i + (j-1)*obj.Nx + (k-1)*obj.Nx*obj.Ny;
                    obj.Neighbours(g).indexes = [g-1, g+1, g+obj.Nx, g-obj.Nx*obj.Ny, g+obj.Nx*obj.Ny];
                end
            end
            j = obj.Ny;
            for k=2:obj.Nz-1
                for i = 2:obj.Nx-1
                    g = i + (j-1)*obj.Nx + (k-1)*obj.Nx*obj.Ny;
                    obj.Neighbours(g).indexes = [g-1, g+1, g-obj.Nx, g-obj.Nx*obj.Ny, g+obj.Nx*obj.Ny];
                end
            end
            k = 1;
            for j=2:obj.Ny-1
                for i=2:obj.Nx-1
                    g = i + (j-1)*obj.Nx + (k-1)*obj.Nx*obj.Ny;
                    obj.Neighbours(g).indexes = [g-1, g+1, g-obj.Nx, g+obj.Nx, g+obj.Nx*obj.Ny];
                end
            end
            k = obj.Nz;
            for j=2:obj.Ny-1
                for i = 2:obj.Nx-1
                    g = i + (j-1)*obj.Nx + (k-1)*obj.Nx*obj.Ny;
                    obj.Neighbours(g).indexes = [g-1, g+1, g-obj.Nx, g+obj.Nx, g-obj.Nx*obj.Ny];
                end
            end
            
            
            % All inner cells (easy ones)
            for k = 2:obj.Nz-1
                for i = 2:obj.Nx-1
                    for j=2:obj.Ny-1
                        g = i + (j-1)*obj.Nx + (k-1)*obj.Nx*obj.Ny;
                        obj.Neighbours(g).indexes = [g-1, g+1, g-obj.Nx, g+obj.Nx, g-obj.Nx*obj.Ny, g+obj.Nx*obj.Ny];
                    end
                end
            end
        end
        function AssignNeighbours2D(obj)
            i = 1;
            for j=1:obj.Ny
                g = i + (j-1)*obj.Nx;
                if j~=1 && j~= obj.Ny
                    obj.Neighbours(g).indexes = [g+1, g-obj.Nx, g+obj.Nx];
                elseif j == 1
                    obj.Neighbours(g).indexes = [g+1, g+obj.Nx];
                elseif j == obj.Ny
                    obj.Neighbours(g).indexes = [g+1, g-obj.Nx];
                end
            end
            i = obj.Nx;
            for j=1:obj.Ny
                g = i + (j-1)*obj.Nx;
                if j~=1 && j~=obj.Ny
                    obj.Neighbours(g).indexes = [g-1, g-obj.Nx, g+obj.Nx];
                elseif j == 1
                    obj.Neighbours(g).indexes = [g-1, g+obj.Nx];
                elseif j == obj.Ny
                    obj.Neighbours(g).indexes = [g-1, g-obj.Nx];
                end
            end
            j = 1;
            for i=2:obj.Nx-1
                g = i + (j-1)*obj.Nx;
                obj.Neighbours(g).indexes = [g-1, g+1, g+obj.Nx];
            end
            j = obj.Ny;
            for i = 2:obj.Nx-1
                g = i + (j-1)*obj.Nx;
                obj.Neighbours(g).indexes = [g-1, g+1, g-obj.Nx];
            end
            
            for i = 2:obj.Nx-1
                for j=2:obj.Ny-1
                    g = i + (j-1)*obj.Nx;
                    obj.Neighbours(g).indexes = [g-1, g+1, g-obj.Nx, g+obj.Nx];
                end
            end
        end
        function AssignNeighbours1D(obj)
            i = 1;
            g = i;
            obj.Neighbours(g).indexes = g+1;
            i = obj.Nx;
            g = i;
            obj.Neighbours(g).indexes = g-1;
            for i = 2:obj.Nx-1
                g = i;
                obj.Neighbours(g).indexes = [g-1, g+1];
            end  
        end
        function AddWells(obj, Inj, Prod)
            for c = 1:obj.N
                Well_Index = [];
                for w = 1:length(Inj)
                    if ~isempty(intersect( Inj(w).Cells , [obj.Children{c,:}] ))
                        Well_Index = [Well_Index; w];
                    end
                end
                for w = length(Inj)+1 : length(Inj)+length(Prod)
                    if ~isempty(intersect( Prod(w-length(Inj)).Cells , [obj.Children{c,:}] ))
                        Well_Index = [Well_Index; w];
                    end
                end
                obj.Wells{c} = Well_Index;
                if isempty(obj.Wells{c})
                    obj.Wells{c} = 0;
                end
            end
        end
    end
end