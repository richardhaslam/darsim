%  Coarse grid class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 26 July 2016
%Last modified: 26 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef coarse_grid < handle
    properties
        N
        Nx
        Ny
        CoarseFactor
        I
        J
        Father
        Active
        Wells
        Neighbours
        Centers   
    end
    methods
        function obj = coarse_grid()
           
        end
        function BuildCoarseGrid(obj, FineGrid)
            %Construct a coarse Grid given coarsening ration and fine Grid
            obj.Nx = FineGrid.Nx/obj.CoarseFactor(1);
            obj.Ny = FineGrid.Ny/obj.CoarseFactor(2);
            obj.N = obj.Nx*obj.Ny;
            %Coordinates of the centres
            obj.I = ones(obj.N, 1);
            obj.J = ones(obj.N, 1);
            Jindexes = ceil(obj.CoarseFactor(2)/2):obj.CoarseFactor(2):FineGrid.Ny;
            for i=1:obj.Ny
                a = obj.Nx*(i-1)+1;
                obj.I(a:a+obj.Nx-1) = ceil(obj.CoarseFactor(1)/2):obj.CoarseFactor(1):FineGrid.Nx;
                obj.J(a:a+obj.Nx-1) = Jindexes(i)*ones(obj.Nx,1);
            end
            obj.Active = zeros(obj.N, 1);
            obj.Wells = zeros(obj.N, 1);
            obj.AssignNeighbours();
        end
        function AssignNeighbours(obj)
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
    end
end