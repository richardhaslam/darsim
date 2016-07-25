%TimeStep Selector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 19 July 2016
%Last modified: 19 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef timestep_selector < handle
    properties
        MinDt
        MaxDt
        ReportDt
        NextDt
        CFL
    end
    methods
        function obj = timestep_selector(cfl)
            obj.MinDt = 10; % in s (hard-coded for now)
            obj.MaxDt = 30*24*3600; % max 30 days
            obj.NextDt = 0;
            obj.CFL = cfl;
        end
        function dt = ComputeStableTimeStep(obj, ProductionSystem, DiscretizationModel, FluidModel, U)
            por = ProductionSystem.Reservoir.por;
            %Void Volume in each cell
            pv = por * DiscretizationModel.ReservoirGrid.Volume;   
            
            %I take the worst possible scenario
            s = Fluid.sr(2):0.01:1-Fluid.sr(1);
            [~, df] =  FluidModel.ComputeFractionalFlow(s);
            dfmax = max(df);
            Uxmax = max(max(abs(U.x)));
            Uymax = max(max(abs(U.y)));
            Lambdax = dfmax * Uxmax;
            Lambday = dfmax * Uymax;
            
            %Compute timestep size
            dtx = obj.CFL*pv/Lambdax;
            dty = obj.CFL*pv/Lambday;
            dt = min(dtx,dty);
        end
        function dt = ChooseTimeStep(obj)
            dt = min([obj.ReportDt, obj.NextDt, obj.MaxDt]);
            dt = max(obj.MinDt, dt);
        end
        function Update(obj, dt, itCount, chops)
            if itCount < 4 && chops < 1
                obj.NextDt = 2*dt;
            elseif itCount > 8 || chops > 1
                obj.NextDt = dt/2;
            else
                obj.NextDt = dt;
            end
        end
    end
end