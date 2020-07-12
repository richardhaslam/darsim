classdef cross_connections < handle
    % Contains the information of non-neighboring connectivities between different media (fracture0matrix and fracture-fracture)
    properties
        Cells
        CI
        T_Geo
        T_Geo_Cond
        UpWind
        U_Geo
    end
    methods
        function obj = cross_connections()
            
        end
        function CrossConnectionsUpWind(obj, P, I)
            [~, n_phases] = size(obj.UpWind);
            for i = 1 : n_phases
                obj.UpWind(:,i) = (P(obj.Cells) - P(I) >= 0);
                obj.U_Geo(:,i) = obj.T_Geo .* (P(obj.Cells) - P(I));
            end
        end
    end
    
end

