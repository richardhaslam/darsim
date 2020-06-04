% Quadratic relative permeability model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 15 July 2016
%Last modified: 15 December 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef relperm_model_table < relperm_model
    properties
        Table
        Type
    end
    methods
        function obj = relperm_model_table(type,table)
            % Each table has 4 columns (Sw, Pc, Krw, Krnw).
            obj.Type = type;
            switch type
                case('Imbibition')
                    obj.Table.Imbibition.TableData = table{1};
                case('Drainage')
                    obj.Table.Draianage.TableData = table{1};
                case('Cyclic')
                    obj.Table.Imbibition.TableData = table{1};
                    obj.Table.Draianage.TableData = table{2};
            end
            obj.SplineCurveFit();
        end
        function SplineCurveFit(obj)
            % Continue here
            Sw = obj.Table.Imbibition(:,1);
            krw = obj.Table.Imbibition(:,3);
            krnw = obj.Table.Imbibition(:,4);
            kr(:,1) = spline(Sw,krw);
            %obj.Table.Imbibition.SplineData = ...
            %kr = coeff(:,1).*Sw.^3 + coeff(:,2).*Sw.^2 + coeff(:,3).*Sw.^1 + coeff(:,4);
        end
        function kr = ComputeRelPerm(obj, Phases, s)
            % Rescale saturations
            S = (s-Phases(1).sr)/(1-Phases(1).sr-Phases(2).sr);
            S = max(S, 0);
            % Phase 1 relative permeability
            kr(:,1) = S.^obj.n(1);
            kr(s < Phases(1).sr, 1) = 0;
            kr(s < Phases(1).sr, 2) = 1;
            % Phase 2 relative permeability
            kr(:,2) = (1-S).^obj.n(2);
            kr(s > 1 - Phases(2).sr, 2) = 0;
            kr(s > 1 - Phases(2).sr, 1) = 1;
        end
        function dkr = ComputeDerivative(obj, Phases, s)
            S = (s-Phases(1).sr)/(1-Phases(1).sr-Phases(2).sr);
            
            dkr(:,1) = (1-Phases(1).sr-Phases(2).sr)^(-1)* obj.n(1) * S.^(obj.n(1)-1);
            dkr(s < Phases(1).sr, 1) = 0;
            dkr(s < Phases(1).sr, 2) = 0;
            dkr(:,2) = -(1-Phases(1).sr-Phases(2).sr)^(-1)*obj.n(2)*(1-S).^(obj.n(2)-1);
            dkr(s > 1 - Phases(2).sr, 2) = 0;
            dkr(s > 1 - Phases(2).sr, 1) = 0;
        end
        function ddkr = ComputeSecondDerivative(obj, Phases, s)
            ddkr(:,1) = ones(length(s), 1) * (1-Phases(1).sr-Phases(2).sr)^(-1)*2;
            ddkr(s < Phases(1).sr, 1) = 0;
            ddkr(s < Phases(1).sr, 2) = 0;
            ddkr(:,2) = ones(length(s), 1) * (1-Phases(1).sr-Phases(2).sr)^(-1)*2;
            ddkr(s > 1 - Phases(2).sr, 2) = 0;
            ddkr(s > 1 - Phases(2).sr, 1) = 0;
        end
    end
end