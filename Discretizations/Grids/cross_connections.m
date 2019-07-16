classdef cross_connections < handle
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    properties
        Cells
        CI
        T_Geo
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

