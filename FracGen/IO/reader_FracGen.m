% Output writer base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 13 July 2016
%Last modified: 13 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef reader_FracGen < handle
    properties
        Directory
        File
        InputMatrix
    end
    methods
        function obj = reader_FracGen(dir, file)
            obj.Directory = dir;
            obj.File = strcat(dir,'/',file);
        end
        function ReadInputFile(obj)
            %ReadInputFile
            fileID = fopen(obj.File, 'r');
            %// Read lines from input file
            matrix = textscan(fileID, '%s', 'Delimiter', '\n');
            obj.InputMatrix = matrix;
            fclose(fileID);
        end
    end
end