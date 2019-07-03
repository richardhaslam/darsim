%fixed TimeStep Selector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Ludovica Delpopolo
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef fixed_timestep_selector < timestep_selector
    properties

    end
    methods

        function dt = ChooseTimeStep(obj)
            %Small time steps for the first steps and then bigger fixed
            %global time steps.
            if obj.Index <= 5
                dt = obj.MinDt;
            else
                dt = obj.MaxDt;

                
            end
%             if (obj.ReportDt == obj.FirstReportDt) && ( dt<obj.PreviousDt || dt<obj.BeforePreviousDt ) && ( obj.NextDt > obj.PreviousDt )
%                 dt = min(obj.ReportDt, obj.MaxDt);
%             end
        end
        function Update(obj, dt, itCount, chops)
            if itCount <= 6 && chops < 1
               obj.NextDt = 2*dt;
            elseif itCount > 12 || chops > 1
               obj.NextDt = dt/2;
            else
               obj.NextDt = dt;
            end
        end
        function UpdateSequential(obj, dt, itOuter ,itTransp, chops)
               obj.NextDt = dt;
        end
    end
end
