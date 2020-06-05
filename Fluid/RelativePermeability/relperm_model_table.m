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
        Sw_from_Table
        krw_from_Table
        krnw_from_Table
        krw_curve
        krnw_curve
    end
    methods
        function obj = relperm_model_table(type,table)
            % Each table has 4 columns (Sw, Pc, Krw, Krnw).
            obj.Type = type;
            switch type
                case('Imbibition')
                    obj.Table.Data = table{1};
                    obj.S_irr(1) = min( obj.Table.Data(:,1) );
                    obj.S_irr(2) = 1 - max( obj.Table.Data(:,1) );
                case('Drainage')
                    obj.Table.Data = table{1};
                    obj.S_irr(2) = min( obj.Table.Data(:,1) );
                    obj.S_irr(1) = 1 - max( obj.Table.Data(:,1) );
                case('Cyclic')
%                     obj.Table.Imbibition.TableData = table{1};
%                     obj.Table.Draianage.TableData = table{2};
                    error('Cyclic rel perm is not implemented yet');
            end
            obj.SplineCurveFit();
        end
        function SplineCurveFit(obj)
            obj.Sw_from_Table = obj.Table.Data(:,1);
            obj.krw_from_Table = obj.Table.Data(:,3);
            obj.krnw_from_Table = obj.Table.Data(:,4);
            % Obtaining the curvature of the rel perm from the tables using spline function
            obj.krw_curve = spline(obj.Sw_from_Table,obj.krw_from_Table);
            obj.krnw_curve = spline(obj.Sw_from_Table,obj.krnw_from_Table);
        end
        function kr = ComputeRelPerm(obj, Phases, s)
            % Rescale saturations
            S = (s-Phases(1).sr)/(1-Phases(1).sr-Phases(2).sr);
            S = max(S, 0);

            % Phase 1 relative permeability
            kr(:,1) = spline(1-obj.Sw_from_Table,obj.krnw_from_Table, S);
            kr(s <= Phases(1).sr, 1) = 0;
            kr(s >= 1-Phases(2).sr, 1) = 1;
            
            % Phase 2 relative permeability
            kr(:,2) = spline(obj.Sw_from_Table,obj.krnw_from_Table, S);
            kr(1-s <= Phases(2).sr, 2) = 0;
            kr(1-s >= 1-Phases(1).sr, 2) = 1;
        end
        function dkr = ComputeDerivative(obj, Phases, s)
            S = (s-Phases(1).sr)/(1-Phases(1).sr-Phases(2).sr);
            
            % Calcuting the first derivative of spline curve using MATLAB built-in function fnder
            dkrw_curve = fnder(obj.krw_curve,1);
            dkrnw_curve = fnder(obj.krnw_curve,1);

            % Obtaining the first derivative of the rel perm using the spline curve derivative
            dkr(:,1) = ppval(dkrw_curve,S);
            dkr(s < Phases(1).sr, 1) = 0;
            dkr(s < Phases(1).sr, 2) = 0;
            
            dkr(:,2) = ppval(dkrnw_curve,S);
            dkr(s > 1 - Phases(2).sr, 2) = 0;
            dkr(s > 1 - Phases(2).sr, 1) = 0;
        end
        function ddkr = ComputeSecondDerivative(obj, Phases, s)
            S = (s-Phases(1).sr)/(1-Phases(1).sr-Phases(2).sr);
            
            % Calcuting the second derivative of spline curve using MATLAB built-in function fnder
            ddkrw_curve = fnder(obj.krw_curve,2);
            ddkrnw_curve = fnder(obj.krnw_curve,2);
            
            % Obtaining the second derivative of the rel perm using the spline curve second derivative
            ddkr(:,1) = ppval(ddkrw_curve,S);
            ddkr(s < Phases(1).sr, 1) = 0;
            ddkr(s < Phases(1).sr, 2) = 0;
            
            ddkr(:,2) = ppval(ddkrnw_curve,S);
            ddkr(s > 1 - Phases(2).sr, 2) = 0;
            ddkr(s > 1 - Phases(2).sr, 1) = 0;
        end
    end
end