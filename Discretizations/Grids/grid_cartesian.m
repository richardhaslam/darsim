%  Cartesian grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 12 July 2016
%Last modified: 14 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef grid_cartesian < grid
    properties
        Nx
        Ny
        Nz
    end
    methods
        function obj = grid_cartesian(nx, ny, nz)
            obj.Nx = nx;
            obj.Ny = ny;
            obj.Nz = nz;
        end
    end
end