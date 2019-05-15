% Buckley-Leverett refined cell selector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef RefCellsSelector_BL < handle
    properties
        tol = 1e-1;
    end
    methods
        function Refinement = SelectRefinedCells(obj, ProductionSystem)
            % We choose to resolve with a fine dt all cells which were
            % flooded during the coarse dt. So we check which cells went
            % from a quasi-zero saturation to something close to the shock
            % one.
            %For now it only works for the reservoir and not for the
            %fractures
            dS = ProductionSystem.Reservoir.State.Properties(obj.key).Value - ...
                        ProductionSystem.Reservoir.State_old.Properties(obj.key).Value;
            Refinement = dS >= obj.tol;
            % P = speye(N);
            % P = P(Cells, :);
        end
    end
end