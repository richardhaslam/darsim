% Productio or Injection curves
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Curve_Inj_Prod < handle
properties
    Phases
    Components
end
methods
    function obj = Curve_Inj_Prod(MaxNTimeSteps, n_phases, n_components, Wells)
        n = length(Wells);
        obj.Phases = zeros(MaxNTimeSteps, n, n_phases);
        obj.Components = zeros(MaxNTimeSteps, n, n_components);
    end
end
end