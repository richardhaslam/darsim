%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Created: 21 March 2016
%Last modified: 15 February 2017 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function creates a vtk comparison file between 2 runs. It can be
%opened in Paraview to create plots of the errors of the different
%variables. 

function VTKErrorFile(Case1, Case2, ComparisonDir, n_files)
%Read the input file to acquire info about geometry etc.
InputFile = strcat(Case1, 'BOHetero.txt');
fileID = fopen(InputFile, 'r');
%// Read lines from input file
inputMatrix = textscan(fileID, '%s', 'Delimiter', '\n');
fclose(fileID);
temp = strfind(inputMatrix{1}, 'TITLE'); % Search a specific string and find all rows containing matches
Problem1 = char(inputMatrix{1}(find(~cellfun('isempty', temp)) + 1));

%Read the input file to acquire info about geometry etc.
InputFile = strcat(Case2, 'BOHetero.txt');
fileID = fopen(InputFile, 'r');
%// Read lines from input file
inputMatrix = textscan(fileID, '%s', 'Delimiter', '\n');
fclose(fileID);
temp = strfind(inputMatrix{1}, 'TITLE'); % Search a specific string and find all rows containing matches
Problem2 = char(inputMatrix{1}(find(~cellfun('isempty', temp)) + 1));

ComparisonFile = strcat(ComparisonDir, Problem2, '_Comparison_');

%%%%%%%%%%%%%PROPERTIEs OF THE RESERVOIR%%%%%%%%%%%%%%%%
temp = strfind(inputMatrix{1}, 'DIMENS'); % Search a specific string and find all rows containing matches
size = find(~cellfun('isempty', temp));
temp = strfind(inputMatrix{1}, 'SPECGRID');
grid = find(~cellfun('isempty', temp));
Grid.Nx = str2double(inputMatrix{1}(grid + 1));
Grid.Ny = str2double(inputMatrix{1}(grid + 2));
Grid.Nz = str2double(inputMatrix{1}(grid + 3));
Grid.Lx = str2double(inputMatrix{1}(size + 1));  %Dimension in x−direction [m]
Grid.Ly = str2double(inputMatrix{1}(size + 2));  %Dimension in y−direction [m]
Grid.h = str2double(inputMatrix{1}(size + 3));  %Reservoir thickness [m]
Grid.dx = Grid.Lx/Grid.Nx;
Grid.dy = Grid.Ly/Grid.Ny;
Grid.dz = Grid.h/Grid.Nz;
Grid.N = Grid.Nx * Grid.Ny * Grid.Nz;
clear temp size grid perm por x

Nstops = n_files;

%% Load data from files
for i=0:Nstops
file = strcat(Case1, 'Output/Solution/',Problem1,'_Sol', num2str(i),'.txt');
Sol_1 = load(file);

file = strcat(Case2, 'Output/Solution/',Problem2,'_Sol', num2str(i),'.txt');
Sol_2 = load(file);

    %Write a VTK file
    fileID = fopen(strcat(ComparisonFile, num2str(i),'.vtk'), 'w');
    fprintf(fileID, '# vtk DataFile Version 2.0\n');
    fprintf(fileID, strcat(Problem2, ' results: Matteo Simulator\n'));
    fprintf(fileID, 'ASCII\n');
    fprintf(fileID, '\n');
    fprintf(fileID, 'DATASET RECTILINEAR_GRID\n');
    fprintf(fileID, 'DIMENSIONS    %d   %d   %d\n', Grid.Nx +1, Grid.Ny+1, 2);
    fprintf(fileID, '\n');
    fprintf(fileID, ['X_COORDINATES ' num2str(Grid.Nx+1) ' float\n']);
    fprintf(fileID, '%f ', 0:Grid.dx:Grid.Lx);
    fprintf(fileID, '\n');
    fprintf(fileID, ['Y_COORDINATES ' num2str(Grid.Ny+1) ' float\n']);
    fprintf(fileID, '%f ', 0:Grid.dy:Grid.Ly);
    fprintf(fileID, '\n');
    fprintf(fileID, 'Z_COORDINATES 2 float\n');
    fprintf(fileID, '%d ', [0 1]);
    fprintf(fileID, '\n');
    fprintf(fileID, '\n');
    fprintf(fileID, 'CELL_DATA   %d\n', Grid.N);
    fprintf(fileID, '\n');
     
    %Pressure error
    PressureError = abs(Sol_2(:,2) - Sol_1(:,2))/(100 - 10);
    PrintScalar2VTK(fileID, PressureError, ' PressureError');
    %Saturation error
    SaturationError = Sol_2(:,3) - Sol_1(:,3);
    PrintScalar2VTK(fileID, SaturationError, ' SaturationError');
    %z error
    SaturationError = Sol_2(:,4) - Sol_1(:,4);
    PrintScalar2VTK(fileID, SaturationError, ' ZError');
end
end
