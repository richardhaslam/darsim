%  Explicit transport solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 26 July 2016
%Last modified: 26 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef explicit_transport_solver < handle
    properties
        itCount
        SystemBuilder
    end
    methods
        function obj = explicit_transport_solver()
            obj.itCount = 1;
        end
        function Solve(obj, ProductionSystem, FluidModel, DiscretizationModel, Formulation, dt)
            disp('explicit transport solver');      
            % 1. Compute fractional flow
            obj.SystemBuilder.ComputePropertiesAndDerivatives(Formulation, ProductionSystem, FluidModel, DiscretizationModel);
            
            % 2. Assemble system matrix
            Formulation.ViscousMatrix(DiscretizationModel.ReservoirGrid);
            
            % 3. Solve explicit transport
            Formulation.UpdateSaturationExplicitly(ProductionSystem, DiscretizationModel, dt);
            
        end
    end
end