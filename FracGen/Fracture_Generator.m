% Fracture Generator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Mousa HosseiniMehr
%TU Delft
%Created: 2017-03-10
%Last modified: 2017-03-10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Fracture_Generator < handle
    properties
        Reader
        Builder
        Simulation
        Writer
    end
    methods
        function obj = Fracture_Generator(Directory, File)
            obj.Reader  = reader_frac(Directory, File);
            obj.Builder = builder_frac();          
        end
        function BuildObjects(obj)
            obj.Builder.FindKeyWords(obj.Reader.InputMatrix);
            obj.Builder.BuildFracInput(obj.Reader.InputMatrix);
            obj.Simulation = obj.Builder.BuildSimulation(obj.Reader.InputMatrix); 
        end
        function PrintInfo(obj)
            disp( ['**************** Generating Fractures for ', obj.Simulation.Domain,' Domain Simulation *****************']);
            disp( char(5));
            disp( 'Reservoir geometry:');
            disp( ['Length: ', num2str(obj.Simulation.Reservoir.LX), ' [m]'] );
            disp( ['Width : ', num2str(obj.Simulation.Reservoir.LY), ' [m]'] );
            if obj.Simulation.Domain == '2D'
                disp( ['Grid  : ', num2str(obj.Simulation.Reservoir.NX), ' x ',  num2str(obj.Simulation.Reservoir.NY)] );
            elseif obj.Simulation.Domain == '3D'
                disp( ['Depth : ', num2str(obj.Simulation.Reservoir.LZ), ' [m]'] );
                disp( ['Grid  : ', num2str(obj.Simulation.Reservoir.NX), ' x ',  num2str(obj.Simulation.Reservoir.NY), ' x ', num2str(obj.Simulation.Reservoir.NZ)] );
            else
                Error( 'The Domain (2D or 3D) is not mentioned correctly in the input file! Is there a typo?' );
            end
            disp( '---------------------------------------------------------' );
        end
        function Run(obj)
            % Calculating the fracture connectivities and populating the fractures with data
            if obj.Simulation.Domain == '2D'
                obj.Simulation.Fractures = FracGen2D(obj);
            elseif obj.Simulation.Domain == '3D'
                obj.Simulation.Fractures = FracGen3D(obj);
            else
                Error( 'The Domain (2D or 3D) is not mentioned correctly in the input file! Is there a typo?' );
            end
            disp( '---------------------------------------------------------' );
            % Write initial state on a file
            % obj.Writer = Plot_VTK('Output/', 'EDFM_3D');
            % obj.Writer.PlotSolution(obj.Simulation);
        end
        function OutputResults(obj)
            Output_File = '../FracGen_IO/Fracture_Output.txt';
            Fracture_Writer(Output_File, obj);
        end
    end
end
