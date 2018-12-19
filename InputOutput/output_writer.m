% Output writer base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef output_writer < handle
    properties
        Directory
        ProblemName
        Plotter
        FormatSol
        FormatInj
        FormatProd
        FormatStats
        FormatTimers
    end
    methods
        function obj = output_writer(dir, problem, n_inj, n_prod, n_timers,  n_stats, ~)
            obj.Directory = strcat(dir, '/Output/');
            obj.ProblemName = problem;
            
            if ~exist(strcat(obj.Directory,'Solution/'), 'dir')
                mkdir(obj.Directory,'Solution/');
            else
                delete(strcat(obj.Directory,'Solution/*.txt'));
            end
            
            %Format for production output
            obj.FormatInj = '%10.5f';
            for i=1:n_inj
                obj.FormatInj = [obj.FormatInj, ' %10.5e'];
                i = i+1;
            end
            obj.FormatInj = [obj.FormatInj, ' %10.5e\n'];
            
            % Format for Injection output
            obj.FormatProd = '%10.5f';
            for i=1:n_prod
                obj.FormatProd = [obj.FormatProd, ' %10.5e'];
                i = i+1;
            end
            obj.FormatProd = [obj.FormatProd, ' %10.5e\n'];
            
            %Format for timers
            obj.FormatTimers = '%10.0f';
            for i = 2:n_timers 
                obj.FormatTimers = [obj.FormatTimers, ' %10.3f'];
            end
            obj.FormatTimers = [obj.FormatTimers, ' %10.3f\n'];
            
            %Format for Stats
            obj.FormatStats = '%10.0f';
            for i = 2:n_stats 
                obj.FormatStats = [obj.FormatStats, ' %10.0f'];
            end
            obj.FormatStats = [obj.FormatStats, ' %10.2f\n'];
            
            %Format for Solution
            obj.FormatSol = '%10.2f';
            for i = 2:2 
                obj.FormatSol = [obj.FormatSol, ' %10.5f'];
            end
            obj.FormatSol = [obj.FormatSol, ' %10.5f\n'];
        end
        function AddPlotter(obj, plotter)
            obj.Plotter = plotter;
        end
        function WriteSolutionOnFile(obj, ProductionSystem, index)
            cells = 1:length(ProductionSystem.Reservoir.State.Properties('P_1').Value);
            fileID = fopen(strcat(obj.Directory, 'Solution/', obj.ProblemName,'_Sol',num2str(index),'.txt'),'w');
            fprintf(fileID, obj.FormatSol, [cells', ProductionSystem.Reservoir.State.Properties('P_1').Value/1e5, ProductionSystem.Reservoir.State.Properties('S_1').Value]');
            for f = 1:ProductionSystem.FracturesNetwork.NumOfFrac
                cells = cells(end)+1:cells(end)+length(ProductionSystem.FracturesNetwork.Fractures(f).State.Properties('P_1').Value);
                fprintf(fileID, obj.FormatSol, [cells', ProductionSystem.FracturesNetwork.Fractures(f).State.Properties('P_1').Value/1e5, ProductionSystem.FracturesNetwork.Fractures(f).State.Properties('S_1').Value]');
            end
            fclose(fileID);
        end
        function WriteWellsData(obj, Time, WellsData, Ndt)
            for i = 1:WellsData.NofPhases
                %Injection
                fileIDI = fopen(strcat(obj.Directory,'Inj_Phase',num2str(i),'.txt'),'w');
                fprintf(fileIDI, obj.FormatInj, [Time(1:Ndt), WellsData.Injection.Phases(2:Ndt+1,:,i), sum(WellsData.Injection.Phases(2:Ndt+1,:,i), 2)]');
                fclose(fileIDI);
                %Production
                fileIDP = fopen(strcat(obj.Directory,'Prod_Phase',num2str(i),'.txt'),'w');
                fprintf(fileIDP, obj.FormatProd, [Time(1:Ndt), WellsData.Production.Phases(2:Ndt+1,:,i), sum(WellsData.Production.Phases(2:Ndt+1,:,i), 2)]');
                fclose(fileIDP);
            end
            for i=1:WellsData.NofComp
                %Injection
                fileIDI = fopen(strcat(obj.Directory,'Inj_Comp',num2str(i),'.txt'),'w');
                fprintf(fileIDI, obj.FormatInj, [Time(1:Ndt), WellsData.Injection.Components(2:Ndt+1,:,i), sum(WellsData.Injection.Components(2:Ndt+1,:,i), 2)]');
                fclose(fileIDI);
                %Production
                fileIDP = fopen(strcat(obj.Directory,'Prod_Comp',num2str(i),'.txt'),'w');
                fprintf(fileIDP, obj.FormatProd, [Time(1:Ndt), WellsData.Production.Components(2:Ndt+1,:,i), sum(WellsData.Production.Components(2:Ndt+1,:,i), 2)]');
                fclose(fileIDP);
            end
        end
        function WriteCouplingStats(obj, CouplingStats, Ndt)
            %Stats
            fileID = fopen(strcat(obj.Directory,'SolverStats.txt'),'w');
            fprintf(fileID, '%10s %10s %10s %10s %10s %10s \n', 'Timestep', 'Chops', 'Iterations', 'CFL Global', 'IterationsLTS', 'CFL Local');
            fprintf(fileID, obj.FormatStats, CouplingStats.StatsMatrix(Ndt));
            fclose(fileID);
            %Timers
            fileID = fopen(strcat(obj.Directory,'Timers.txt'),'w');
            fprintf(fileID, '%10s %10s %10s %10s %10s\n', '# Timestep', 'Total Time', 'Construct', 'Solve', 'Flash');
            fprintf(fileID, obj.FormatTimers, CouplingStats.TimersMatrix(Ndt));
            fclose(fileID);
        end
    end
end