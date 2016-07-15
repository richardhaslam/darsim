%  Cartesian grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 12 July 2016
%Last modified: 12 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef cartesian_grid < grid
    properties
        Nx
        Ny
        Nz
    end
    methods
        function obj = cartesian_grid(nx, ny)
            obj.Nx = nx;
            obj.Ny = ny;
            obj.Nz = 1;
            obj.N = nx * ny * 1;
        end
    end
end