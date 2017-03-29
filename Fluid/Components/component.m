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
        MM % Molecular mass
        kval0 % Ref kvalue for correlations
    end
    methods
        function AddCompProperties(obj, mm, kval)
           obj.MM = mm;
           obj.kval0 = kval;
        end
    end
end