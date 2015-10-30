%OUTPUT Statistics on a file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (strcmp(Strategy, 'Sequential')==1)
    Statistics = [Sequential.ImplicitSolver.timestep; Sequential.ImplicitSolver.Chops; Sequential.ImplicitSolver.Newtons];
    fileID = fopen(strcat(Directory,'SeqStat.txt'),'w');
    fileID2 = fopen(strcat(Directory,'SeqTimings.txt'),'w');
    fprintf(fileID,'%6s %12s %12s\n','Timestep','# Chops', '# Newtons');
    fprintf(fileID,'%6.0f %12.0f %12.0f\n', Statistics);
    fclose(fileID);
    Timesteps=linspace(1,Ndt-1,Ndt-1);
    Timers = [Timesteps; TimerTimestep(1:Ndt-1)'; TimerPressure(1:Ndt-1)'; TimerBalance(1:Ndt-1)'; TimerSaturation(1:Ndt-1)'];
    fprintf(fileID2,'%6s %12s %12s %12s  %12s\n','Timestep', 'Total Time','Pressure-Solver', 'Balance-check', 'Transport-Solver');
    fprintf(fileID2,'%6.0f %12.3f %12.3f %12.3f %12.3f\n', Timers);
    fclose(fileID2);
else
    Statistics = [FIM.timestep(1:Ndt-1), FIM.Chops(1:Ndt-1), FIM.Iter(1:Ndt-1)];
    fileID = fopen(strcat(Directory,'FIMStat.txt'),'w');
    fileID2 = fopen(strcat(Directory,'FIMTimings.txt'),'w');
    fileID3 = fopen(strcat(Directory,'FIMOilProd.txt'),'w');
    fileID4 = fopen(strcat(Directory,'FIMSaturation.txt'),'w');
    fileID5 = fopen(strcat(Directory,'FIMPressure.txt'),'w');
    %fprintf(fileID,'%6s %12s %12s\n','Timestep', '# Chops', '# Iterations');
    fprintf(fileID,'%6.0f %12.0f %12.0f\n', Statistics');
    fclose(fileID);
    Timesteps = linspace(1,Ndt-1,Ndt-1);
    Timers = [Timesteps; TimerTimestep(1:Ndt-1)'; TimerConstruct(1:Ndt-1)';  TimerSolve(1:Ndt-1)'];
    fprintf(fileID2,'%6s %12s %12s %12s\n','Timestep', 'Total Time', 'Construct', 'Solve');
    fprintf(fileID2,'%6.0f %12.3f %12.3f %12.3f\n', Timers);
    fclose(fileID2);
    fprintf(fileID3,'%12.3f %12.3f\n', [CumulativeTime(1:Ndt)'; Prod.oil(1:Ndt)']);
    fclose(fileID3);
    fprintf(fileID4,'%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n', Saturations');
    fclose(fileID4);
    fprintf(fileID5,'%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n', Pressures');
    fclose(fileID5);
end