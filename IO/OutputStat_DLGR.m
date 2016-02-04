%DLGR: OUTPUT Statistics on a file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (strcmp(Strategy, 'Sequential')==1)
    disp('No DLGR for Sequential Strategy');
else
    Statistics = [FIM.timestep(1:Ndt-1), FIM.Chops(1:Ndt-1), FIM.Iter(1:Ndt-1), FIM.ActiveCells(1:Ndt-1, :)];
    fileID1 = fopen(strcat(Directory,'DLGRStat.txt'),'w');
    fileID2 = fopen(strcat(Directory,'DLGRTimings.txt'),'w');
    %fileID3 = fopen(strcat(Directory,'DLGROilProd.txt'),'w');
    fileID4 = fopen(strcat(Directory,'Saturation.txt'),'w');
    fileID5 = fopen(strcat(Directory,'Pressure.txt'),'w');
    %fprintf(fileID,'%6s %12s %12s %12.s\n','Timestep', '# Chops', '# Iterations', '# Active Cells');
    fprintf(fileID1,'%6.0f %12.0f %12.0f %12.0f %12.0f %12.0f\n', Statistics');
    fclose(fileID1);
    Timers = [FIM.timestep(1:Ndt-1)'; TimerTimestep(1:Ndt-1)'; TimerRP(1:Ndt-1)'; TimerConstruct(1:Ndt-1)';  TimerSolve(1:Ndt-1)'];
    fprintf(fileID2,'%6s %12s %12s %12s %12s\n','Timestep', 'Total Time', 'R and P', 'Jacobian', 'Solve');
    fprintf(fileID2,'%6.0f %12.3f %12.3f %12.3f %12.3f\n', Timers);
    fclose(fileID2);
    %fprintf(fileID3,'%12.3f %12.3f\n', [CumulativeTime(1:Ndt)'; Prod.oil(1:Ndt)']);
    %fclose(fileID3);
    fprintf(fileID4,'%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n', Saturations');
    fclose(fileID4);
    fprintf(fileID5,'%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n', Pressures');
    fclose(fileID5);
end