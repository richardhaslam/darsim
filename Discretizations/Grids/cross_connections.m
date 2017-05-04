classdef cross_connections < handle
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    properties
        Cells
        T_Geo
        UpWind
        U_Geo
    end
    
    methods
        function CrossConnectionsUpWind(obj, ProductionSystem, DiscretizationModel, If1_Local, Formulation)
            obj.UpWind = false(length(obj.Cells), Formulation.NofPhases);
            obj.U_Geo = zeros(length(obj.Cells), Formulation.NofPhases);
            for i = 1 : Formulation.NofPhases
                If1_Global = If1_Local + DiscretizationModel.ReservoirGrid.N;
                Index_frac1_Local = DiscretizationModel.Index_Global_to_Local(If1_Global);
                Pf = ProductionSystem.FracturesNetwork.Fractures(Index_frac1_Local.f).State.Properties(['P_',num2str(i)]).Value(Index_frac1_Local.g);
                % UpWind of frac-matrix
                indices_m = obj.Cells( obj.Cells <= DiscretizationModel.ReservoirGrid.N );
                Pm = ProductionSystem.Reservoir.State.Properties(['P_',num2str(i)]).Value(indices_m);    
                obj.UpWind(1:length(indices_m),i) = Pm - Pf >= 0;
                Obj.U_Geo(1:length(indices_m),i) = obj.T_Geo(1:length(indices_m)) .* (Pm - Pf);
                % Upwind of frac-farc
                indices_f = obj.Cells( obj.Cells > DiscretizationModel.ReservoirGrid.N );
                if ~isempty(indices_f)
                    for n = 1:length(indices_f)
                        If2_Global = indices_f(n); % Global indices of the other fractures' cells if any
                        If2_Local = If2_Global - DiscretizationModel.ReservoirGrid.N;
                        Index_frac2_Local = DiscretizationModel.Index_Global_to_Local(If2_Global);
                        Pf_other = ProductionSystem.FracturesNetwork.Fractures(Index_frac2_Local.f).State.Properties(['P_',num2str(i)]).Value(Index_frac2_Local.g);
                        obj.UpWind(length(indices_m)+n,i) = Pf_other - Pf >= 0;
                        obj.U_Geo(length(indices_m)+n,i) = obj.T_Geo(length(indices_m)+n) * (Pf_other - Pf);
                    end
                end
            end
        end
    end
    
end

