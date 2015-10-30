%READ INPUT File
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%ReadInputFile
fileID = fopen(InputFile, 'r');
%// Read lines from input file
C = textscan(fileID, '%s', 'Delimiter', '\n');
fclose(fileID);

%Reservoir Properties
%// Search a specific string and find all rows containing matches
grid = strfind(C{1}, 'GRID');
rows = find(~cellfun('isempty', C));

%Fluid Properties

%Simulator Settings
