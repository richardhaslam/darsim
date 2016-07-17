% Wells 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 13 July 2016
%Last modified: 13 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef wells < handle
    properties
        Inj
        Prod
        NofInj
        NofProd
    end
    methods
        function AddInjector(obj, Injector)
            obj.Inj = [obj.Inj, Injector];
        end
        function AddProducer(obj, Producer)
            obj.Prod = [obj.Prod, Producer];
        end
    end
end