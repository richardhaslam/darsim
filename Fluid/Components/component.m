% Component class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 15 July 2016
%Last modified: 13 April 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef component < matlab.mixin.Heterogeneous & handle
    properties
        MM % Molecular mass
        kval0 % Ref kvalue for correlations
        Pcrit
        Tcrit
        w % acentric factor 
    end
    methods
        function AddCompProperties(obj, Prop)
           obj.MM = Prop(1);
           obj.kval0 = Prop(2);
           if length(Prop) > 2
               obj.Tcrit = Prop(3);
               obj.Pcrit = Prop(4);
               obj.w = Prop(5);
           end
        end
    end
end