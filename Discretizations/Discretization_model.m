%  Discretization model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 12 July 2016
%Last modified: 12 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Discretization_model < handle
    properties
        ReservoirGrid
        FracturesGrid
    end
    methods
        function obj = Discretization_model(nx, ny, nz)
            obj.ReservoirGrid = cartesian_grid(nx, ny, nz);
        end
        function DefinePerforatedCells(obj, Wells)
            % Has to be improved for Diagonal wells (maybe using trajectories)
            % Injectors
            for i = 1:Wells.NofInj
                I = Wells.Inj(i).Coord(1,1):1:Wells.Inj(i).Coord(1,2);
                Y = Wells.Inj(i).Coord(2,1):1:Wells.Inj(i).Coord(2,2);
                Wells.Inj(i).Cells = I + (Y-1)*obj.ReservoirGrid.Nx;
            end
            
            % Producers
            for i = 1:Wells.NofProd
                I = Wells.Prod(i).Coord(1,1):1:Wells.Prod(i).Coord(1,2);
                Y = Wells.Prod(i).Coord(2,1):1:Wells.Prod(i).Coord(2,2);
                Wells.Prod(i).Cells = I + (Y-1)*obj.ReservoirGrid.Nx;
            end
        end
        function AverageMassOnCoarseBlocks(obj, Status, FluidModel, Formulation)
            % virtual call
        end
    end
end