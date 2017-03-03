% Status 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Matteo Cusini
%TU Delft
%Created: 14 July 2016
%Last modified: 3 March 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef status < matlab.mixin.Copyable
    % Defines the State of the Reservoir
properties
    Properties = containers.Map % contains all existing properties
    T
    p
    pc
    S
    x
    z
    rho
    rhoT
    ni
end
methods
    function Initialize(obj, N, nf, ncomp)
        % Create objects
        obj.p = ones(N, 1);
        obj.pc = zeros(N, 1);
        obj.z = zeros(N, ncomp);
        obj.S = zeros(N, 1);
        obj.ni = 0.5*ones(N, 1);
        obj.x = zeros(N, nf*ncomp);
        obj.x(:,1) = 1;
        obj.x(:,end) = 1;
        obj.rho = zeros(N, nf);
        obj.rhoT = zeros(N, 1); 
    end
    function AddProperties(obj, FluidModel, N)
        % Add properties to property map
        %% Pressure
        Prop = property(N, 1);
        obj.Properties('Pressure') = Prop;
        %% Density
        obj.Properties('rho_1') = Prop;
        obj.Properties('rhoT') = Prop;
        switch (FluidModel.name)
            case ('SinglePhase')
               % No more properties to be added
            case ('Immiscible')
                %% Saturation and Pc
                for i=2:FluidModel.NofPhases
                    obj.Properties(['S_', num2str(i)]) = Prop;
                    obj.Properties(['rho_', num2str(i)]) = Prop;
                end
                obj.Properties('Pc') = Prop;
            otherwise
                %% Saturation and Pc
                for i=1:FluidModel.NofPhases
                    obj.Properties(['S_', num2str(i)]) = Prop;
                    obj.Properties(['rho_', num2str(i)]) = Prop;
                    obj.Properties(['ni_', num2str(i)]) = Prop;
                end
                obj.Properties('Pc') = Prop;
                %% Compositional stuff
                for i=1:FluidModel.NofComp
                    obj.Properties(['z_', num2str(i)]) = Prop;
                end
                Prop = property(N, FluidModel.NofPhases * FluidModel.NofComp);
                obj.Properties('x') = Prop;
                
        end
    end
    function AssignInitialValues(obj, VarNames, VarValues)
        for i=1:length(VarNames)
            Prop = obj.Properties(char(VarNames{i}));
            Prop.Value = ones(length(Prop.Value), 1) * VarValues(i);
        end
    end
    function UpdatePressure(obj)
    end
    function UpdateSaturation(obj)
    end
    function UpdateZ(obj)
    end
end
end