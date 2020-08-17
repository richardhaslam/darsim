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
    end
    methods
        function obj = reader(dir, file)
            obj.Directory = dir;
            obj.File = strcat(dir,'/',file);
        end
    end
    methods (Abstract)
        obj = ReadInputFile(obj);
    end
end