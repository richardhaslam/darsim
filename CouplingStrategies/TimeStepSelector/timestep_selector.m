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
            rho = zeros(1, FluidModel.NofPhases);
            for i=1:FluidModel.NofPhases
                rho(i) = max(ProductionSystem.Reservoir.State.Properties(['rho_', num2str(i)]).Value);
                Mob(:, i) = rho(i) * Mob(:, i);
                dMob(:, i) = rho(i) * dMob(:, i); 
            end
            num = rho(:, 1)  .* Mob(:, 1);
            dnum = rho(:, 1) .* dMob(:,1);
            den = sum(rho .* Mob, 2);
            dden = sum(rho .* dMob, 2);  
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
            dtx = obj.CFL*pv*max(rho)/Lambdax;
            dty = obj.CFL*pv*max(rho)/Lambday;
            dtz = obj.CFL*pv*max(rho)/Lambdaz;
            dt = min([dtx, dty, dtz, obj.ReportDt, obj.MaxDt]);
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
    end
end