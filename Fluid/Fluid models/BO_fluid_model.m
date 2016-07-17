% Immiscible Fluid model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 14 July 2016
%Last modified: 15 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef BO_fluid_model < fluid_model
    properties
        Pref % Pref for Rs computation
    end
    methods
        function obj = BO_fluid_model(n_phases, n_comp)
            obj@fluid_model(n_phases, n_comp);
        end
        function Status = InitializeReservoir(obj, Status)
        end
        function Inj = InitializeInjectors(obj, Inj)
        end
    end
end