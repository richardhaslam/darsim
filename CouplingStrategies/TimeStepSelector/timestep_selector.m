%TimeStep Selector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef timestep_selector < handle
    properties
		TotalTime		 
        MinDt
        MaxDt
        ReportDt
		FirstReportDt			 
        NextDt
		PreviousDt = 0;
        BeforePreviousDt = 0;			   
        CFL
        Index
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
            dMob = FluidModel.ComputeDMobDS(s);
           
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
															
            dtx = obj.CFL*pv/Lambdax;
            dty = obj.CFL*pv/Lambday;
            dtz = obj.CFL*pv/Lambdaz;
            % This means that the CFL should be the max CFL you want to use
            dt = min([dtx, dty, dtz, obj.ReportDt, obj.MaxDt, obj.NextDt]);
        end
        function CreateCFLfile(obj, ProductionSystem, DiscretizationModel, df, U, dt)
            por = ProductionSystem.Reservoir.Por;
            %Void Volume in each cell
            pv = por * DiscretizationModel.ReservoirGrid.Volume;   
            
            Ux = abs(U.x(2:end,:,:) + U.x(1:end-1,:,:)/2);
            Uy = abs(U.y(:,2:end,:) + U.y(:,1:end-1,:)/2);
            Uz = abs(U.z(:,:,2:end) + U.z(:,:,1:end-1)/2);

            Lambdax = abs(df) .* Ux(:);
            Lambday = abs(df) .* Uy(:);
            Lambdaz = abs(df) .* Uz(:);
            
            cfl = dt * (Lambdax + Lambday + Lambdaz) / pv;
           
        end
        function dt = ChooseTimeStep(obj)
            if obj.ReportDt <= 0
                obj.ReportDt = obj.MaxDt;
            end

            dt = min([obj.ReportDt, obj.NextDt, obj.MaxDt]);

            %dt = min([obj.ReportDt, obj.NextDt, obj.MaxDt]);
            % If the previous "obj.ReportDt" forces the previous "dt" to be
            % smaller than its two previous consecutive timesteps
            % ("obj.PreviousDt" and "obj.BeforePreviousDt"), and if
            % timestep has not been chopped due to convergence issues, the
            % current  timestep is too small, and it can be as big as
            % "obj.ReportDt".
            
            if (obj.ReportDt == obj.FirstReportDt) && ( dt<obj.PreviousDt || dt<obj.BeforePreviousDt ) && ( obj.NextDt > obj.PreviousDt )
                dt = min(obj.ReportDt, obj.MaxDt);
            end
        end
        function Update(obj, dt, itCount, chops)
            if itCount <= 10 && chops < 1
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
        function Time = Sec2DHMS(obj, TimeInSec)
            Time.Days  = floor( TimeInSec / (24*3600) );
            TimeInSec  = TimeInSec - Time.Days * 24*3600;
            
            Time.Hours = floor( TimeInSec / 3600 );
            TimeInSec  = TimeInSec - Time.Hours * 3600;
            
            Time.Minutes = floor( TimeInSec / 60 );
            
            Time.Seconds  = TimeInSec - Time.Minutes * 60;
        end
    end
end
