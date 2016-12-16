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
        function Initialize(obj, ProductionSystem, Formulation)
            obj.ReservoirGrid.Initialize(ProductionSystem.Reservoir);
            % Perforated cells
            obj.DefinePerforatedCells(ProductionSystem.Wells);
            % . Assign Depth
            obj.ReservoirGrid.ComputeDepth(Formulation.GravityModel.alpha);
        end
        function DefinePerforatedCells(obj, Wells)
            % Has to be improved for Diagonal wells (maybe using trajectories)
            % Injectors
            for i = 1:Wells.NofInj
                I = Wells.Inj(i).Coord(1,1):1:Wells.Inj(i).Coord(1,2);
                J = Wells.Inj(i).Coord(2,1):1:Wells.Inj(i).Coord(2,2);
                K = Wells.Inj(i).Coord(3,1):1:Wells.Inj(i).Coord(3,2);
                Wells.Inj(i).Cells = I + (J-1)*obj.ReservoirGrid.Nx + (K-1)*obj.ReservoirGrid.Nx*obj.ReservoirGrid.Ny;
                Wells.Inj(i).ResizeObjects(length(Wells.Inj(i).Cells));
            end
            
            % Producers
            for i = 1:Wells.NofProd
                I = Wells.Prod(i).Coord(1,1):1:Wells.Prod(i).Coord(1,2);
                J = Wells.Prod(i).Coord(2,1):1:Wells.Prod(i).Coord(2,2);
                K = Wells.Prod(i).Coord(3,1):1:Wells.Prod(i).Coord(3,2);
                Wells.Prod(i).Cells = I + (J-1)*obj.ReservoirGrid.Nx + (K-1)*obj.ReservoirGrid.Nx*obj.ReservoirGrid.Ny;
                Wells.Prod(i).ResizeObjects(length(Wells.Prod(i).Cells));
            end
        end
        function AverageMassOnCoarseBlocks(obj, Status, FluidModel, Formulation)
            % virtual call
        end
    end
end