% Reader base class
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
        function ReadInputFile(obj, Builder)
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
            fractured = find(~cellfun('isempty', temp));
            if isempty(fractured)
                fractured = 0;
            else
                fractured = 1;
            end
            if fractured 
                FractureFile = strcat(obj.Directory, '/Fracture_Output.txt');
                fileID = fopen(FractureFile, 'r');
				fprintf('Reading fracture file ...');
                obj.FractureMatrix = textscan(fileID, '%s', 'Delimiter', '\n');
                fclose(fileID);
				fprintf('--> Successful.\n\n');
            end
        end
    end
end