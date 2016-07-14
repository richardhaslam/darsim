% Output txt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 4 July 2016
%Last modified: 11 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef output_writer_txt < output_writer
    properties
        FormatInj
        FormatProd
        FormatStats
        FormatTimers
    end
    methods
        function obj = output_writer_txt(dir, problem, n_inj, n_prod, n_timers,  n_stats)
            obj@output_writer(dir, problem);
            %Format for production output
            obj.FormatInj = '%10.2f';
            for i=1:n_inj
                obj.FormatInj = [obj.FormatInj, ' %10.2e'];
                i = i+1;
            end
            obj.FormatInj = [obj.FormatInj, ' %10.2e\n'];
            
            % Format for Injection output
            obj.FormatProd = '%10.2f';
            for i=1:n_prod
                obj.FormatProd = [obj.FormatProd, ' %10.2e'];
                i = i+1;
            end
            obj.FormatProd = [obj.FormatProd, ' %10.2e\n'];
            
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
            obj.FormatStats = [obj.FormatStats, ' %10.0f\n'];
        end
        function WriteSummary(obj, Summary)
            obj.WriteWellsData(Summary.Time, Summary.WellsData, Summary.NumberTimeSteps);
            obj.WriteCouplingStats(Summary.CouplingStats, Summary.NumberTimeSteps);
        end
        function WriteWellsData(obj, Time, WellsData, Ndt)
            fileIDI(1) = fopen(strcat(obj.Directory,'Inj_Phase1.txt'),'w');
            fileIDI(2) = fopen(strcat(obj.Directory,'Inj_Phase2.txt'),'w');
            fileIDI(3) = fopen(strcat(obj.Directory,'Inj_Comp1.txt'),'w');
            fileIDI(4) = fopen(strcat(obj.Directory,'Inj_Comp2.txt'),'w');
            fileIDP(1) = fopen(strcat(obj.Directory,'Prod_Phase1.txt'),'w');
            fileIDP(2) = fopen(strcat(obj.Directory,'Prod_Phase2.txt'),'w');
            fileIDP(3) = fopen(strcat(obj.Directory,'Prod_Comp1.txt'),'w');
            fileIDP(4) = fopen(strcat(obj.Directory,'Prod_Comp2.txt'),'w');
            for i = 1:2
                %Injection
                fprintf(fileIDI(i), obj.FormatInj, [Time(1:Ndt), WellsData.Injection.Phases(2:Ndt+1,:,i), sum(WellsData.Injection.Phases(2:Ndt+1,:,i), 2)]');
                fclose(fileIDI(i));
                fprintf(fileIDI(i+2), obj.FormatInj, [Time(1:Ndt), WellsData.Injection.Components(2:Ndt+1,:,i), sum(WellsData.Injection.Components(2:Ndt+1,:,i), 2)]');
                fclose(fileIDI(i+2));
                %Production
                fprintf(fileIDP(i), obj.FormatProd, [Time(1:Ndt), WellsData.Production.Phases(2:Ndt+1,:,i), sum(WellsData.Production.Phases(2:Ndt+1,:,i), 2)]');
                fclose(fileIDP(i));
                fprintf(fileIDP(i+2), obj.FormatProd, [Time(1:Ndt), WellsData.Production.Components(2:Ndt+1,:,i), sum(WellsData.Production.Components(2:Ndt+1,:,i), 2)]');
                fclose(fileIDP(i+2));
            end
        end
        function WriteCouplingStats(obj, CouplingStats, Ndt)
            %Stats
            fileID = fopen(strcat(obj.Directory,'SolverStats.txt'),'w');
            fprintf(fileID, obj.FormatStats, CouplingStats.StatsMatrix(Ndt));
            fclose(fileID);
            %Timers
            fileID = fopen(strcat(obj.Directory,'Timers.txt'),'w');
            fprintf(fileID, obj.FormatTimers, CouplingStats.TimersMatrix(Ndt));
            fclose(fileID);
        end
    end
end