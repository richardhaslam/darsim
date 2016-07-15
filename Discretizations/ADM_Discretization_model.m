%  ADM discretization model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 12 July 2016
%Last modified: 12 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef ADM_Discretization_model < Discretization_model
    properties
        ADMSettings
        CoarseGrid
        ADMGrid
    end
    methods
        function obj = ADM_Discretization_model(Nx, Ny, ADMSettings)
            obj@Discretization_model(Nx, Ny);
            obj.ADMSettings = ADMSettigns;
    end
end