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
        dx
        dy
        dz
        Ax
        Ay
        Volume
    end
    methods
        function obj = grid_cartesian(nx, ny, nz)
            obj.Nx = nx;
            obj.Ny = ny;
            obj.Nz = nz;
            obj.N = nx * ny * nz;
        end
        function Initialize(obj, Reservoir)
            obj.dx = Reservoir.Length/obj.Nx;
            obj.dy = Reservoir.Width/obj.Ny;
            obj.dz = Reservoir.Depth/obj.Nz;
            obj.Ax = obj.dy * obj.dz;
            obj.Ay = obj.dx * obj.dz;
            obj.Volume = obj.dx * obj.dy * obj.dz;
        end
    end
end