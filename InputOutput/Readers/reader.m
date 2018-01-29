% Reader base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 13 July 2016
%Last modified: 14 December 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef reader < handle
    properties
        Directory
        File
        PermDirectory
    end
    methods
        function obj = reader(dir, file, permdirectory)
            obj.Directory = dir;
            obj.File = strcat(dir,'/',file);
            obj.PermDirectory = permdirectory;
        end
    end
    methods (Abstract)
        obj = ReadInputFile(obj);
    end
end