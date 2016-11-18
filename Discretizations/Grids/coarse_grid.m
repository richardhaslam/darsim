%  Coarse grid class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 26 July 2016
%Last modified: 26 July 2016
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
    end
    methods
        function obj = coarse_grid()
           
        end
        function BuildCoarseGrid(obj, FineGrid)
            % Construct a coarse Grid given coarsening ration and fine Grid
            obj.Nx = FineGrid.Nx/obj.CoarseFactor(1);
            obj.Ny = FineGrid.Ny/obj.CoarseFactor(2);
            obj.Nz = FineGrid.Nz/obj.CoarseFactor(3);
            obj.N = obj.Nx*obj.Ny*obj.Nz;
            % Coordinates of the centres
            obj.I = ones(obj.N, 2);
            obj.J = ones(obj.N, 2);
            obj.K = ones(obj.N, 2);
            Jindexesf = ceil(obj.CoarseFactor(2)/2):obj.CoarseFactor(2):FineGrid.Ny;
            Jindexes = 1:1:obj.Ny;
            for k = 1:obj.Nz
                for i=1:obj.Ny
                    a = obj.Nx*(i-1) + (k-1)*obj.Nx*obj.Ny + 1;
                    obj.I(a:a+obj.Nx-1, 2) = ceil(obj.CoarseFactor(1)/2):obj.CoarseFactor(1):FineGrid.Nx;
                    obj.J(a:a+obj.Nx-1, 2) = Jindexesf(i)*ones(obj.Nx,1);
                    obj.I(a:a+obj.Nx-1, 1) = 1:1:obj.Nx;
                    obj.J(a:a+obj.Nx-1, 1) = Jindexes(i)*ones(obj.Nx, 1);
                end
            end
            Kindexes = ceil(obj.CoarseFactor(3)/2):obj.CoarseFactor(3):FineGrid.Nz;
            for i=1:obj.Nz
                obj.K((i-1)*obj.Nx*obj.Ny+1:i*obj.Nx*obj.Ny, 2) = Kindexes(i) * ones(obj.Nx*obj.Ny, 1);
                obj.K((i-1)*obj.Nx*obj.Ny+1:i*obj.Nx*obj.Ny, 1) = i*ones(obj.Nx*obj.Ny, 1);
            end
            obj.Active = zeros(obj.N, 1);
            obj.Wells = zeros(obj.N, 1);
            obj.Fathers = zeros(obj.N, 1);
            obj.AssignNeighbours();
        end
        function AssignNeighbours(obj)
            % Let s do the 8 corners separetely
            % 1.
            i = 1;
            j = 1;
            k = 1;
            g = i + (j-1)*obj.Nx + (k-1)*obj.Nx*obj.Ny;
            obj.Neighbours(g).indexes = [g+1, g+obj.Nx, g+obj.Nx*obj.Ny];
            % 2.
            i = 1;
            j = 1;
            k = obj.Nz;
            g = i + (j-1)*obj.Nx + (k-1)*obj.Nx*obj.Ny;
            obj.Neighbours(g).indexes = [g+1, g+obj.Nx, g-obj.Nx*obj.Ny];
            % 3.
            i = 1;
            j = obj.Ny;
            k = 1;
            g = i + (j-1)*obj.Nx + (k-1)*obj.Nx*obj.Ny;
            obj.Neighbours(g).indexes = [g+1, g-obj.Nx, g+obj.Nx*obj.Ny];
            % 4.
            i = 1;
            j = obj.Ny;
            k = obj.Nz;
            g = i + (j-1)*obj.Nx + (k-1)*obj.Nx*obj.Ny;
            obj.Neighbours(g).indexes = [g+1, g-obj.Nx, g-obj.Nx*obj.Ny];
            % 5.
            i = obj.Nx;
            j = 1;
            k = 1;
            g = i + (j-1)*obj.Nx + (k-1)*obj.Nx*obj.Ny;
            obj.Neighbours(g).indexes = [g-1, g+obj.Nx, g+obj.Nx*obj.Ny];
            % 6.
            i = obj.Nx;
            j = 1;
            k = obj.Nz;
            g = i + (j-1)*obj.Nx + (k-1)*obj.Nx*obj.Ny;
            obj.Neighbours(g).indexes = [g-1, g+obj.Nx, g-obj.Nx*obj.Ny];
            % 7.
            i = obj.Nx;
            j = obj.Ny;
            k = 1;
            g = i + (j-1)*obj.Nx + (k-1)*obj.Nx*obj.Ny;
            obj.Neighbours(g).indexes = [g-1, g-obj.Nx, g+obj.Nx*obj.Ny];
            % 8.
            i = obj.Nx;
            j = obj.Ny;
            k = obj.Nz;
            g = i + (j-1)*obj.Nx + (k-1)*obj.Nx*obj.Ny;
            obj.Neighbours(g).indexes = [g-1, g-obj.Nx, g-obj.Nx*obj.Ny];
            
            % Let's do the borders first
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
    end
end