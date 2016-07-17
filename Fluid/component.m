% Component class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 15 July 2016
%Last modified: 15 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef component < handle
    properties
        Tb % Reference density
        b  % 
    end
    methods
        function obj = component(Tb, b)
           obj.Tb = Tb;
           obj.b = b;
        end
       
    end
end