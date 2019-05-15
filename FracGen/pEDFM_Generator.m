% Fracture Generator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Mousa HosseiniMehr
%TU Delft
%Created: 2017-03-10
%Last modified: 2017-03-10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef pEDFM_Generator < handle
    properties
        Reader
        Builder
        Geometry
        Writer 
    end
    methods
        function obj = pEDFM_Generator(Directory, File)
            obj.Reader  = reader_FracGen(Directory, File);
            obj.Builder = builder_FracGen();
            obj.Writer = writer_FracGen(Directory,'Fracture_Output.txt');
        end
        function BuildObjects(obj)
            obj.Builder.FindKeyWords(obj.Reader.InputMatrix);
            obj.Builder.BuildFracInput(obj.Reader.InputMatrix);
            obj.Geometry = obj.Builder.BuildGeometry(obj.Reader.InputMatrix); 
        end
        function PrintInfo(obj)
            obj.Geometry.ReservoirGrid.PrintInfo(obj.Geometry.Domain);
            obj.Geometry.FracturesGrid.PrintInfo();
        end
        function Run(obj)
            obj.Geometry.PrintInfo();
            obj.Geometry.Compute_EDFM_Connectivities();
            if strcmp(obj.Geometry.Type,'pEDFM')
                obj.Geometry.Add_pEDFM_Connectivities();
            end
            close all;
            obj.Geometry.PlotFractures();
        end
        function OutputResults(obj)
            obj.Writer.VTK_Plotter = plotVTK_FracGen(strcat(obj.Writer.Directory,'/Output'), 'pEDFM_3D');
            obj.Writer.VTK_Plotter.PlotSolution(obj.Geometry);
            obj.Writer.PrintInfo();
            obj.Writer.WriteOutputTextFile(obj.Geometry);
        end
    end
end
