% Output writer base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 13 July 2016
%Last modified: 13 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef reader < handle
    properties
        Directory
        File
        InputMatrix
        SettingsMatrix
        FractureMatrix
    end
    methods
        function obj = reader(dir, file)
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
            %ReadSettingFile
            SettingsFile = strcat(obj.Directory, '/SimulatorSettings.txt');
            fileID = fopen(SettingsFile, 'r');
            matrix = textscan(fileID, '%s', 'Delimiter', '\n');
            obj.SettingsMatrix = matrix;
            fclose(fileID);
            %ReadSettingFile
            temp = strfind(obj.InputMatrix{1}, 'FRACTURED');
            if str2double(  obj.InputMatrix{1}{ find(~cellfun('isempty', temp))+1 }  ) == 1
                FractureFile = strcat(obj.Directory, '/Fracture_Output.txt');
                fileID = fopen(FractureFile, 'r');
                obj.FractureMatrix = textscan(fileID, '%s', 'Delimiter', '\n');
                fclose(fileID);
            end
        end
    end
end