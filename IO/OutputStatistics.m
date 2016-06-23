%OUTPUT Statistics on a file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Created: 2015
%Last modified: 18 May 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[row, columns] = size(Pressures);

%Format for solution output
format = '%10.5f';
for i=1:columns-2
    format = [format, ' %10.5f'];
    i = i+1;
end
format = [format, ' %10.5f\n'];

%Format for production output
format2 = '%10.5f';
for i=1:length(Prod) - 1
    format2 = [format2, ' %10.5f'];
    i = i+1;
end
format2 = [format2, ' %10.5f\n'];

if (strcmp(Strategy, 'Sequential')==1)
    Statistics = [Sequential.ImplicitSolver.timestep; Sequential.ImplicitSolver.Chops; Sequential.ImplicitSolver.Newtons];
    fileID = fopen(strcat(Directory,'SeqStat.txt'),'w');
    fileID2 = fopen(strcat(Directory,'SeqTimings.txt'),'w');
    fileID3a = fopen(strcat(Directory,'SeqNwProd.txt'),'w');
        fileID3b = fopen(strcat(Directory,'SeqWProd.txt'),'w');
    fileID4 = fopen(strcat(Directory,'SeqSaturation.txt'),'w');
    fileID5 = fopen(strcat(Directory,'SeqPressure.txt'),'w');
    fprintf(fileID,'%6s %12s %12s\n','Timestep','# Chops', '# Newtons');
    fprintf(fileID,'%6.0f %12.0f %12.0f\n', Statistics);
    fclose(fileID);
    Timesteps=linspace(1,Ndt-1,Ndt-1);
    Timers = [Timesteps; TimerTimestep(1:Ndt-1)'; TimerPressure(1:Ndt-1)'; TimerBalance(1:Ndt-1)'; TimerSaturation(1:Ndt-1)'];
    fprintf(fileID2,'%6s %12s %12s %12s  %12s\n','Timestep', 'Total Time','Pressure-Solver', 'Balance-check', 'Transport-Solver');
    fprintf(fileID2,'%6.0f %12.3f %12.3f %12.3f %12.3f\n', Timers);
    fclose(fileID2);
else
    if (ADMSettings.active == 0)
        Statistics = [FIM.timestep(1:Ndt-1), FIM.Chops(1:Ndt-1), FIM.Iter(1:Ndt-1)];
        fileID = fopen(strcat(Directory,'FIMStat.txt'),'w');
        fileID2 = fopen(strcat(Directory,'FIMTimings.txt'),'w');
        fileID3a = fopen(strcat(Directory,'FIMNwProd.txt'),'w');
        fileID3b = fopen(strcat(Directory,'FIMWProd.txt'),'w');
        fileID4 = fopen(strcat(Directory,'FIMSaturation.txt'),'w');
        fileID5 = fopen(strcat(Directory,'FIMPressure.txt'),'w');
        %fprintf(fileID,'%6s %12s %12s\n','Timestep', '# Chops', '# Iterations');
        fprintf(fileID,'%6.0f %12.0f %12.0f\n', Statistics');
        fclose(fileID);
        Timesteps = linspace(1,Ndt-1,Ndt-1);
        Timers = [Timesteps; TimerTimestep(1:Ndt-1)'; TimerConstruct(1:Ndt-1)';  TimerSolve(1:Ndt-1)'; TimerInner(1:Ndt-1)'];
        fprintf(fileID2,'%6s %12s %12s %12s %12s\n','Timestep', 'Total Time', 'Construct', 'Solve', 'Flash');
        fprintf(fileID2,'%6.0f %12.3f %12.3f %12.3f %12.3f\n', Timers);
        fclose(fileID2);
        %fprintf(fileID3,'%12.3f %12.3f\n', [CumulativeTime(1:Ndt)'; Prod.oil(1:Ndt)']);
        %fclose(fileID3);
    else
        Statistics = [FIM.timestep(1:Ndt-1), FIM.Chops(1:Ndt-1), FIM.Iter(1:Ndt-1), FIM.ActiveCells(1:Ndt-1, :)];
        fileID1 = fopen(strcat(Directory,'ADMStat.txt'),'w');
        fileID2 = fopen(strcat(Directory,'ADMTimings.txt'),'w');
        fileID3a = fopen(strcat(Directory,'ADMNwProd.txt'),'w');
        fileID3b = fopen(strcat(Directory,'ADMWProd.txt'),'w');
        fileID4 = fopen(strcat(Directory,'ADMSaturation.txt'),'w');
        fileID5 = fopen(strcat(Directory,'ADMPressure.txt'),'w');
        %fprintf(fileID,'%6s %12s %12s %12.s\n','Timestep', '# Chops', '# Iterations', '# Active Cells');
        fprintf(fileID1,'%6.0f %12.0f %12.0f %12.0f %12.0f %12.0f\n', Statistics');
        fclose(fileID1);
        Timers = [FIM.timestep(1:Ndt-1)'; TimerTimestep(1:Ndt-1)'; TimerRP(1:Ndt-1)'; TimerConstruct(1:Ndt-1)';  TimerSolve(1:Ndt-1)'; TimerInner(1:Ndt-1)'];
        fprintf(fileID2,'%6s %12s %12s %12s %12s %12s\n','Timestep', 'Total Time', 'R and P', 'Jacobian', 'Solve', 'Flash');
        fprintf(fileID2,'%6.0f %12.3f %12.3f %12.3f %12.3f %12.3f\n', Timers);
        fclose(fileID2);
        
    end
end
fprintf(fileID3a, format2, NwProduction(:,1:Ndt));
fclose(fileID3a);
fprintf(fileID3b, format2, WProduction(:,1:Ndt));
fclose(fileID3b);
fprintf(fileID4, format, Saturations');
fclose(fileID4);
fprintf(fileID5,format, Pressures');
fclose(fileID5);