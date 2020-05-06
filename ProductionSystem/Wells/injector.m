% Injector 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef injector < handle
    properties
        Coord
        Cells
        p % injection pressure
        qv % total volumetric rate
        T % injection temperature
        QPhases
        QComponents
        Qh % enthalpy flux
        PI
        BHPDepth
        z
        S
        ni
        x
        x2
        rho
        Mob
        h % enthalphy of injection
        BC_Formulation
        
    end
    methods
        function obj = injector(PI, coord, n_phases)
            obj.PI = PI;
            obj.Coord = coord;
            obj.rho = ones(1, n_phases);
            obj.qv = 0;
        end
        function ResizeObjects(obj, n)
            obj.p =  ones(n,1) * obj.p;
            n_phases = length(obj.rho);
            obj.rho = ones(n, n_phases);
            obj.qv = ones(n, 1) * obj.qv/n;
        end
        function residual = InjectorResidual(obj, residual)
            residual = sum(obj.QPhases);
        end
    end
end