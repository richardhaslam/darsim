% Gravity model class for CornerPointGrid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: 
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef gravity_model_CornerPointGrid < gravity_model
    properties
    end
    methods
        function initialize_rhoint(obj, Grid, n_phases, f)
            for i=1:n_phases
                obj.RhoInt{i, f+1} = zeros(size(Grid.Trans));
            end
        end
        function ComputeInterfaceDensities(obj, Grid, Status, f)
            [n_phases, ~] = size(obj.RhoInt);
            for i=1:n_phases
                rho_g = obj.g * Status.Properties(['rho_', num2str(i)]).Value;
                rho_g_1 = rho_g( Grid.CornerPointGridData.Internal_Face.CellNeighbor1Index );
                rho_g_2 = rho_g( Grid.CornerPointGridData.Internal_Face.CellNeighbor2Index );
                obj.RhoInt{i, f+1} = (rho_g_1 + rho_g_2) / 2;
            end
        end
    end
end

