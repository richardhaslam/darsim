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
            
            FormatSol = '%10.2f\n';
            fileID = fopen(strcat('../Input/SPE10T1/Output/','Solution/','CFL',num2str(obj.Index),'.txt'),'w');
            fprintf(fileID, FormatSol, cfl);
        end
        function dt = ChooseTimeStep(obj)
            if obj.ReportDt <= 0
                obj.ReportDt = obj.MaxDt;
            end
            %dt = obj.NextDt;
            if obj.Index <= 5
                dt = obj.MinDt;
                %dt = min([obj.ReportDt, obj.NextDt, obj.MaxDt]);
            else
                dt = obj.MaxDt;
            end
            %dt = max(obj.MinDt, dt);
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
            % We check for number of outer iterations and nubmer of
            % iterations of last transport solve.
%             if itOuter <= 2 && itTransp < 6 && chops < 1
%                obj.NextDt = 2*dt;
%             elseif itTransp > 10 || chops > 1
%                obj.NextDt = dt/2;
%             else
               obj.NextDt = dt;
%             end
        end
    end
end