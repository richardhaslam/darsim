% Reader for eclipse input files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef reader_eclipse < handle
    properties
        Directory
        File
        PermDirectory
    end
    methods
        function obj = reader_eclipse(dir, file)
            obj.Directory = dir;
            obj.File = strcat(dir,'/',file);
        end
    end
end