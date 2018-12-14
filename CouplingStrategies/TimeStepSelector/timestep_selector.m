%TimeStep Selector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
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
        function obj = timestep_selector(cfl, mindt, maxdt)
            obj.MinDt  = mindt*24*3600; 
            obj.MaxDt  = maxdt*24*3600; 
            obj.NextDt = mindt*24*3600;
            obj.CFL    = cfl;
        end
        function dt = StableTimeStep(obj, ProductionSystem, DiscretizationModel, FluidModel, U)
            por = ProductionSystem.Reservoir.Por;
            %Void Volume in each cell
            pv = por * DiscretizationModel.ReservoirGrid.Volume;   
            
            % It's an approxmiation of the real CFL condition
            %I take the worst possible scenario
            s = [FluidModel.Phases(1).sr:0.01:1-FluidModel.Phases(1).sr]';
            Mob = FluidModel.ComputePhaseMobilities(s);
            dMob = FluidModel.DMobDS(s);
           
            num =  Mob(:, 1);
            dnum = dMob(:,1);
            den = sum(Mob, 2);
            dden = sum(dMob, 2);  
            df = (dnum .* den - dden .* num) ./ den.^2;
            dfmax = max(df);
            Uxmax = max(max(max(abs(U.x))));
            Uymax = max(max(max(abs(U.y))));
            Uzmax = max(max(max(abs(U.z))));
            Lambdax = dfmax * Uxmax;
            Lambday = dfmax * Uymax;
            Lambdaz = dfmax * Uzmax;
            
            %Compute timestep size
            % I multiply by mean(rho) coz U is the mass flux
            dtx = obj.CFL*pv/Lambdax;
            dty = obj.CFL*pv/Lambday;
            dtz = obj.CFL*pv/Lambdaz;
            % This means that the CFL should be the max CFL you want to use
            dt = min([dtx, dty, dtz, obj.ReportDt, obj.MaxDt, obj.NextDt]);
        end
        function dt = ChooseTimeStep(obj)
            if obj.ReportDt <= 0
                obj.ReportDt = obj.MaxDt;
            end
            dt = min([obj.ReportDt, obj.NextDt, obj.MaxDt]);
            %dt = max(obj.MinDt, dt);
        end
        function Update(obj, dt, itCount, chops)
            if itCount <= 4 && chops < 1
               obj.NextDt = 2*dt;
            elseif itCount > 12 || chops > 1
               obj.NextDt = dt/2;
            else
               obj.NextDt = dt;
            end
        end
        function UpdateSequential(obj, dt, itOuter ,itTransp, chops)
            % We check for number of outer iterations and nubmer of
            % iterations of last transport solve.
            if itOuter <= 2 && itTransp < 6 && chops < 1
               obj.NextDt = 2*dt;
            elseif itTransp > 10 || chops > 1
               obj.NextDt = dt/2;
            else
               obj.NextDt = dt;
            end
        end
    end
end