% Reader base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
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
            
            if ~strcmp('/',permdirectory) && ~strcmp('\',permdirectory)
                permdirectory = strcat(permdirectory,'/');
            end
            obj.PermDirectory = permdirectory;
        end
    end
    methods (Abstract)
        obj = ReadInputFile(obj);
    end
end