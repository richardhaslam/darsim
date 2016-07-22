% Wells data base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 4 July 2016
%Last modified: 8 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef wells_data < handle 
    properties
        NofPhases
        NofComp
        Injection
        Production
    end
    methods
        function obj = wells_data(MaxNTimeSteps, n_phases, n_components, Wells)
            obj.Injection = Curve_Inj_Prod(MaxNTimeSteps, n_phases, n_components, Wells.Inj);
            obj.Production = Curve_Inj_Prod(MaxNTimeSteps, n_phases, n_components, Wells.Prod);
            obj.NofPhases = n_phases;
            obj.NofComp = n_components;
        end
        function UpdateInjectionCurve(obj, Ndt, Inj, dT)
            for w = 1:length(Inj)
                for i=1:obj.NofPhases
                    obj.Injection.Phases(Ndt, w, i) = obj.Injection.Phases(Ndt-1, w, i) + sum(Inj(w).QPhases(:,i), 1)*dT;
                end
                for i=1:obj.NofComp
                    obj.Injection.Components(Ndt, w, i) = obj.Injection.Components(Ndt-1, w, i) + sum(Inj(w).QComponents(:,i), 1)*dT;
                end
            end
        end
        function UpdateProductionCurve(obj, Ndt, Prod, dT)
            for w = 1:length(Prod)
                for i=1:obj.NofPhases
                    obj.Production.Phases(Ndt, w, i) = obj.Production.Phases(Ndt-1, w, i) - sum(Prod(w).QPhases(:,i), 1)*dT;
                end
                for i=1:obj.NofComp
                    obj.Production.Components(Ndt, w, i) = obj.Production.Components(Ndt-1, w, i) - sum(Prod(w).QComponents(:, i), 1)*dT;
                end
            end
        end
    end
end