%TimeStep Selector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 19 July 2016
%Last modified: 26 July 2016
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
            obj.MinDt = 100; % in s (hard-coded for now)
            obj.MaxDt = 30*24*3600; % max 30 days
            obj.NextDt = 0.01*24*3600;
            obj.CFL = cfl;
        end
        function dt = StableTimeStep(obj, ProductionSystem, DiscretizationModel, FluidModel, U)
            por = ProductionSystem.Reservoir.Por;
            %Void Volume in each cell
            pv = por * DiscretizationModel.ReservoirGrid.Volume;   
            
            %I take the worst possible scenario
            s = [FluidModel.Phases(1).sr:0.01:1-FluidModel.Phases(1).sr]';
            Mob = FluidModel.ComputePhaseMobilities(s);
            dMob = FluidModel.DMobDS(s);
            df = (dMob(:,1) .* sum(Mob, 2) - sum(dMob, 2) .* Mob(:,1)) ./ sum(Mob, 2).^2;
            dfmax = max(df);
            Uxmax = max(max(max(abs(U.x))));
            Uymax = max(max(max(abs(U.y))));
            Uzmax = max(max(max(abs(U.z))));
            Lambdax = dfmax * Uxmax;
            Lambday = dfmax * Uymax;
            Lambdaz = dfmax * Uzmax;
            
            %Compute timestep size
            dtx = obj.CFL*pv/Lambdax;
            dty = obj.CFL*pv/Lambday;
            dtz = obj.CFL*pv/Lambdaz;
            dt = min([dtx, dty, dtz, obj.ReportDt, obj.MaxDt]);
        end
        function dt = ChooseTimeStep(obj)
            dt = min([obj.ReportDt, obj.NextDt, obj.MaxDt]);
            %dt = max(obj.MinDt, dt);
        end
        function Update(obj, dt, itCount, chops)
            if itCount <= 7 && chops < 1
                obj.NextDt = 2*dt;
            elseif itCount > 12 || chops > 1
                obj.NextDt = dt/2;
            else
                obj.NextDt = dt;
            end
        end
    end
end