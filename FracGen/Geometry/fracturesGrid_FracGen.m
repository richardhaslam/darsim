% class fracturers for DARSim2FracGen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Mousa HosseiniMehr
%TU Delft
%Created: 2016-07-12
%Last modified: 2017-03-01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef fracturesGrid_FracGen < handle
    properties
        N_Frac
        N_total
        Fracture
        Epsilon
        % The level of accuracy for distance comparison between the points is not
        % zero, but a value close to zero (Epsilon). This is needed because of
        % machine's non-zero accuracy
    end
    methods
        function obj = fracturesGrid_FracGen(epsilon)
            %%
            obj.N_Frac = 0;
            obj.N_total = 0;
            obj.Epsilon = epsilon;
            obj.Fracture = fracture_FracGen();
        end
        function BuildFracturesGrid(obj, ReservoirGrid, Frac_Input)
            % Initializing total number of fracture grids
            obj.N_total = 0;
            %% Fractures Construction with Geometry Input
            f = 0;
            for F = 1 : size(Frac_Input,2)
                % "F" is counter for all the fractures in the input file (whether
                % active or not), but "f" is the counter for all the active fractures.
                if Frac_Input(F).isActive
                    f = f + 1;
                    
                    % Building the fracture plane
                    obj.Fracture(f) = fracture_FracGen();
                    obj.Fracture(f).BuildFracture(Frac_Input(F), ReservoirGrid, obj.Epsilon, f);
                    
                    % Cumulative number of all fractures cells
                    obj.N_total = obj.N_total + obj.Fracture(f).N_Length_AB*obj.Fracture(f).N_Width_AD;
                end
            end
            obj.N_Frac = length(obj.Fracture);
            obj.Fracture = obj.Fracture';
            if f == 0
                error('No fracture exists. At least one fracture should be active. Check the input file!');
            end
        end
        function PrintInfo(obj)
            %%
            for f = 1 : obj.N_Frac
                obj.Fracture(f).PrintInfo(f);
            end
            fprintf('\n');
            fprintf('Total number of fractures: %3.0f\n',obj.N_Frac);
            fprintf('Total number of fractures cells: %3.0f\n',obj.N_total);
            fprintf('---------------------------------------------------------\n');
        end
    end
end