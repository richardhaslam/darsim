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
format_s = '%10.4e';
for i=1:columns-2
    format_s = [format_s, ' %10.4e'];
    i = i+1;
end
format_s = [format_s, ' %10.4e\n'];

format_p = '%10.5e';
for i=1:columns-2
    format_p = [format_p, ' %10.5e'];
    i = i+1;
end
format_p = [format_p, ' %10.5e\n'];

%Format for production output
format2 = '%10.2f';
for i=1:length(Prod)
    format2 = [format2, ' %10.2e'];
    i = i+1;
end
format2 = [format2, ' %10.2e\n'];

%Format for production output
format3 = '%10.2f';
for i=1:length(Inj)
    format3 = [format3, ' %10.2e'];
    i = i+1;
end
format3 = [format3, ' %10.2e\n'];

if (strcmp(Strategy, 'Sequential')==1)
    Statistics = [Sequential.ImplicitSolver.timestep; Sequential.ImplicitSolver.Chops; Sequential.ImplicitSolver.Newtons];
    fileID = fopen(strcat(Directory,'SeqStat.txt'),'w');
    fileID2 = fopen(strcat(Directory,'SeqTimings.txt'),'w');
    fileID3a = fopen(strcat(Directory,'Seq_W_Prod.txt'),'w');
    fileID3b = fopen(strcat(Directory,'Seq_Nw_Prod.txt'),'w');
    fileID3c = fopen(strcat(Directory,'Seq_Comp1_Prod.txt'),'w');
    fileID3d = fopen(strcat(Directory,'Seq_Comp2_Prod.txt'),'w');
    fileID4a = fopen(strcat(Directory,'Seq_W_Inj.txt'),'w');
    fileID4b = fopen(strcat(Directory,'Seq_Nw_Inj.txt'),'w');
    fileID4c = fopen(strcat(Directory,'Seq_Comp1_Inj.txt'),'w');
    fileID4d = fopen(strcat(Directory,'Seq_Comp2_Inj.txt'),'w');
    fileID5 = fopen(strcat(Directory,'SeqSaturation.txt'),'w');
    fileID6 = fopen(strcat(Directory,'SeqPressure.txt'),'w');
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
        fileID3a = fopen(strcat(Directory,'FIM_W_Prod.txt'),'w');
        fileID3b = fopen(strcat(Directory,'FIM_Nw_Prod.txt'),'w');
        fileID3c = fopen(strcat(Directory,'FIM_Comp1_Prod.txt'),'w');
        fileID3d = fopen(strcat(Directory,'FIM_Comp2_Prod.txt'),'w');
        fileID4a = fopen(strcat(Directory,'FIM_W_Inj.txt'),'w');
        fileID4b = fopen(strcat(Directory,'FIM_Nw_Inj.txt'),'w');
        fileID4c = fopen(strcat(Directory,'FIM_Comp1_Inj.txt'),'w');
        fileID4d = fopen(strcat(Directory,'FIM_Comp2_Inj.txt'),'w');
        fileID5 = fopen(strcat(Directory,'FIMSaturation.txt'),'w');
        fileID6 = fopen(strcat(Directory,'FIMPressure.txt'),'w');
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
        fileID3a = fopen(strcat(Directory,'ADM_W_Prod.txt'),'w');
        fileID3b = fopen(strcat(Directory,'ADM_Nw_Prod.txt'),'w');
        fileID3c = fopen(strcat(Directory,'ADM_Comp1_Prod.txt'),'w');
        fileID3d = fopen(strcat(Directory,'ADM_Comp2_Prod.txt'),'w');
        fileID4a = fopen(strcat(Directory,'ADM_W_Inj.txt'),'w');
        fileID4b = fopen(strcat(Directory,'ADM_Nw_Inj.txt'),'w');
        fileID4c = fopen(strcat(Directory,'ADM_Comp1_Inj.txt'),'w');
        fileID4d = fopen(strcat(Directory,'ADM_Comp2_Inj.txt'),'w');
        fileID5 = fopen(strcat(Directory,'ADMSaturation.txt'),'w');
        fileID6 = fopen(strcat(Directory,'ADMPressure.txt'),'w');
        %fprintf(fileID,'%6s %12s %12s %12.s\n','Timestep', '# Chops', '# Iterations', '# Active Cells');
        fprintf(fileID1,'%6.0f %12.0f %12.0f %12.0f %12.0f %12.0f\n', Statistics');
        fclose(fileID1);
        Timers = [FIM.timestep(1:Ndt-1)'; TimerTimestep(1:Ndt-1)'; TimerRP(1:Ndt-1)'; TimerConstruct(1:Ndt-1)';  TimerSolve(1:Ndt-1)'; TimerInner(1:Ndt-1)'];
        fprintf(fileID2,'%6s %12s %12s %12s %12s %12s\n','Timestep', 'Total Time', 'R and P', 'Jacobian', 'Solve', 'Flash');
        fprintf(fileID2,'%6.0f %12.3f %12.3f %12.3f %12.3f %12.3f\n', Timers);
        fclose(fileID2);
        
    end
end

%Production
fprintf(fileID3a, format2, [Production.time(1:Ndt); Production.Phase.W(:,1:Ndt); sum(Production.Phase.W(:,1:Ndt), 1)]);
fclose(fileID3a);
fprintf(fileID3b, format2, [Production.time(1:Ndt); Production.Phase.Nw(:,1:Ndt); sum(Production.Phase.Nw(:,1:Ndt), 1)]);
fclose(fileID3b);
fprintf(fileID3c, format2, [Production.time(1:Ndt); Production.Component.z1(:,1:Ndt); sum(Production.Component.z1(:,1:Ndt), 1)]);
fclose(fileID3c);
fprintf(fileID3d, format2, [Production.time(1:Ndt); Production.Component.z2(:,1:Ndt); sum(Production.Component.z2(:,1:Ndt), 1)]);
fclose(fileID3d);


%Injection
fprintf(fileID4a, format3, [Injection.time(1:Ndt); Injection.Phase.W(:,1:Ndt); sum(Injection.Phase.W(:,1:Ndt), 1)]);
fclose(fileID4a);
fprintf(fileID4b, format3, [Injection.time(1:Ndt); Injection.Phase.Nw(:,1:Ndt); sum(Injection.Phase.Nw(:,1:Ndt), 1)]);
fclose(fileID4b);
fprintf(fileID4c, format3, [Injection.time(1:Ndt); Injection.Component.z1(:,1:Ndt); sum(Injection.Component.z1(:,1:Ndt), 1)]);
fclose(fileID4c);
fprintf(fileID4d, format3, [Injection.time(1:Ndt); Injection.Component.z2(:,1:Ndt); sum(Injection.Component.z2(:,1:Ndt),1)]);
fclose(fileID4d);

fprintf(fileID5, format_s, Saturations');
fclose(fileID5);
fprintf(fileID6,format_p, Pressures');
fclose(fileID6);