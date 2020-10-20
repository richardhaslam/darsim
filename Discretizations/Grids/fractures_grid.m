%  fractures grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Mousa HosseiniMehr and Matteo Cusini
%TU Delft
%Created: 21 March 2017
%Last modified: 21 March 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef fractures_grid < handle
properties
    Nfrac
    N
    Grids
end
methods
    function obj = fractures_grid(n_frac)
        obj.Nfrac = n_frac;
        obj.N = zeros(n_frac, 1);
        obj.Grids = cartesian_grid.empty;
    end
    function AddGrid(obj, Grid, i)
        obj.Grids(i,1) = Grid;
        obj.N(i,1) = Grid.N;
    end
end
end