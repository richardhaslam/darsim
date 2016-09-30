% Component class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 15 July 2016
%Last modified: 28 September 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef component < matlab.mixin.Heterogeneous & handle
    properties
        Tb % Bubble point temperature
        b  % 
        MM % Molecular mass
    end
    methods
        function AddCompProperties(obj, Tb, b, mm)
           obj.Tb = Tb;
           obj.b = b;
           obj.MM = mm;
        end
    end
end