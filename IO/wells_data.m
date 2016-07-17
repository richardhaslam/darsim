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
        Injection
        Production
    end
    methods
        function obj = wells_data(MaxNTimeSteps, n_phases, n_components, Wells)
            obj.Injection = Curve_Inj_Prod(MaxNTimeSteps, n_phases, n_components, Wells.Inj);
            obj.Production = Curve_Inj_Prod(MaxNTimeSteps, n_phases, n_components, Wells.Prod);
        end
        function UpdateInjectionCurve(obj, Ndt, Inj, dT)
            for w = 1:length(Inj)
                obj.Injection.Phases(Ndt, w, 1) = obj.Injection.Phases(Ndt-1, w, 1) + Inj(w).qw*dT;
                obj.Injection.Phases(Ndt, w, 2) = obj.Injection.Phases(Ndt-1, w, 2) + Inj(w).qnw*dT;
                obj.Injection.Components(Ndt, w, 1) = obj.Injection.Components(Ndt-1, w, 1) + Inj(w).qz1*dT;
                obj.Injection.Components(Ndt, w, 2) = obj.Injection.Components(Ndt-1, w, 2) + Inj(w).qz2*dT;
            end
        end
        function UpdateProductionCurve(obj, Ndt, Prod, dT)
            for w = 1:length(Prod)
                obj.Production.Phases(Ndt, w, 1) = obj.Production.Phases(Ndt-1, w, 1) - Prod(w).qw*dT;
                obj.Production.Phases(Ndt, w, 2) = obj.Production.Phases(Ndt-1, w, 2) - Prod(w).qnw*dT;
                obj.Production.Components(Ndt, w, 1) = obj.Production.Components(Ndt-1, w, 1) - Prod(w).qz1*dT;
                obj.Production.Components(Ndt, w, 2) = obj.Production.Components(Ndt-1, w, 2) - Prod(w).qz2*dT;
            end
        end
    end
end