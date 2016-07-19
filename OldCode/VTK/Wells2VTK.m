%Output wells in VTK format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Created: 27 May 2016
%Last modified: 27 May 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Wells2VTK(Grid, Inj, Prod, Directory, Problem)
%InJectors
for i=1:length(Inj)
    name = strcat(Directory,'/VTK/',Problem,num2str(i),'Inj','.vtk');
    WriteAWell(Inj(i), name, Grid);
end

%Producers
for i=1:length(Prod)
    name = strcat(Directory,'/VTK/',Problem,num2str(i),'Prod','.vtk');
    WriteAWell(Prod(i), name, Grid);
end

end

function WriteAWell(Well, name, Grid)
fileID = fopen(name, 'w');
fprintf(fileID, '# vtk DataFile Version 2.0\n');
fprintf(fileID, 'Matteo Simulator: Well s file\n');
fprintf(fileID, '\n');
fprintf(fileID, 'ASCII\n');
fprintf(fileID, '\n');
fprintf(fileID, 'DATASET POLYDATA\n');
fprintf(fileID, ['POINTS ', num2str(length(Well.cells) + 1), ' float\n']);
%loop over perforations
z1 = Grid.h/2;
z2 = 2*Grid.h;
for c=1:length(Well.cells)
    dummy = mod(Well.cells(c), Grid.Nx);
    if dummy == 0
        i = Grid.Nx;
    else
        i = dummy;
    end
    j = (Well.cells(c) - i)/Grid.Nx + 1;
    x = (i - 1)*Grid.dx + Grid.dx/2;
    y = (j - 1)*Grid.dy + Grid.dy/2;
    fprintf(fileID, '%10.5f %10.5f %10.5f\n', x, y, z1);
end
fprintf(fileID, '%10.5f %10.5f %10.5f\n', x, y, z2);
fprintf(fileID, ['LINES ', num2str(length(Well.cells)), ' ',num2str(3*(length(Well.cells))), '\n']);
for c = 1:length(Well.cells)
    fprintf(fileID, '%d %d %d\n', [2, c-1, c]);
end
fclose(fileID);
end