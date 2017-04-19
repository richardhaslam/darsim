% Status 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Matteo Cusini
%TU Delft
%Created: 14 July 2016
%Last modified: 9 March 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef status < handle
    % Defines the State of the reservoir or of the fractures
properties
    Properties
    T
end
methods
    function obj = status()
        obj.Properties = containers.Map;
    end
    function AddProperties(obj, FluidModel, N)
        % Add properties to property map
        % Total density
        obj.Properties('rhoT') = property(N, 1, false, 0, 2000);
        switch (FluidModel.name)
            case ('SinglePhase')
               % No more properties to be added
               for i=1:FluidModel.NofPhases
                    obj.Properties(['P_', num2str(i)]) = property(N, 1, true, 1e7, 2e7);
                    obj.Properties(['S_', num2str(i)]) = property(N, 1, false, 0, 1);
                    obj.Properties(['rho_', num2str(i)]) = property(N, 1, false, 0, 2000);
                end
            case ('Immiscible')
                %% Saturation and Pc
                for i=1:FluidModel.NofPhases
                    obj.Properties(['P_', num2str(i)]) = property(N, 1, true, 1e7, 2e7);
                    obj.Properties(['S_', num2str(i)]) = property(N, 1, false, 0, 1);
                    obj.Properties(['rho_', num2str(i)]) = property(N, 1, false, 0, 2000);
                end
                obj.Properties('Pc') = property(N, 1, false, 1e3, 1e6);
            otherwise
                %% Saturation and Pc
                for i=1:FluidModel.NofPhases
                    obj.Properties(['P_', num2str(i)]) = property(N, 1, true, 1e7, 2e7);
                    obj.Properties(['S_', num2str(i)]) = property(N, 1, true, 0, 1);
                    obj.Properties(['rho_', num2str(i)]) = property(N, 1, false, 0, 2000);
                    obj.Properties(['ni_', num2str(i)]) = property(N, 1, false, 0, 1);
                end
                obj.Properties('Pc') = property(N, 1, false, 1e3, 1e6);
                %% Compositional stuff
                for i=1:FluidModel.NofComp
                    obj.Properties(['z_', num2str(i)]) = property(N, 1, true, 0, 1);
                    for j=1:FluidModel.NofPhases
                        obj.Properties(['x_', num2str(i), 'ph', num2str(j)]) = property(N, 1, false, 0, 1, i==j);
                    end
                end
        end
    end
    function AssignInitialValues(obj, VarNames, VarValues)
        for i=1:length(VarNames)
            Prop = obj.Properties(char(VarNames{i}));
            Prop.Value = ones(length(Prop.Value), 1) * VarValues(i);
        end
    end
    function CopyProperties(obj, Source)
        Names = Source.Properties.keys;
        N_prop = double(Source.Properties.Count);
        for i = 1:N_prop
            obj.Properties([Names{i}]) = property(1, 1);
            temp = obj.Properties(Names{i});
            temp.Plot = Source.Properties(Names{i}).Plot;
            temp.Value = Source.Properties(Names{i}).Value;
            temp.Valmax = Source.Properties(Names{i}).Valmax;
            temp.Valmin = Source.Properties(Names{i}).Valmin; 
        end
    end
end
end