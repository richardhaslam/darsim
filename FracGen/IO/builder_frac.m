% Builder Builds all objects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 13 July 2016
%Last modified: 4 March 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef builder_frac < handle
    properties
        dim
        size
        grid
        compareAccuracy
        fracStart
        fracEnd
        NrOfAllFracs
        FracInput
        plotting
    end
    methods
        function FindKeyWords(obj, inputMatrix)
            
            %%%%%%%%%%%%%PROPERTIES OF THE RESERVOIR%%%%%%%%%%%%%%%%
            % Search a specific string and find all rows containing matches
            temp = strfind(inputMatrix{1}, 'DOMAIN');
            obj.dim = find(~cellfun('isempty', temp));
            temp = strfind(inputMatrix{1}, 'DIMENS');
            obj.size = find(~cellfun('isempty', temp));
            temp = strfind(inputMatrix{1}, 'SPECGRID');
            obj.grid = find(~cellfun('isempty', temp));
            
            temp = strfind(inputMatrix{1}, 'CompareAccuracy');
            obj.compareAccuracy = find(~cellfun('isempty', temp));
            
            %%%%%%%%%%%%%PROPERTIES OF THE Fractures%%%%%%%%%%%%%%%%
            temp = strfind(inputMatrix{1}, 'FRAC_START');
            obj.fracStart = find(~cellfun('isempty', temp));
            
            temp = strfind(inputMatrix{1}, 'FRAC_END');
            obj.fracEnd = find(~cellfun('isempty', temp));
            obj.NrOfAllFracs = obj.fracEnd - obj.fracStart - 1;
            
        end
        
        function BuildFracInput(obj, inputMatrix)
            obj.FracInput = struct();
            for f = 1 : obj.NrOfAllFracs
                frac_temp = strsplit(inputMatrix{1}{obj.fracStart+f},' , ');
                
                %%% Check if fracture is used or not %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                obj.FracInput(f).isActive       = str2double( frac_temp{01} );
                %%% Central Coordinate of the fracture plate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                frac_temp{02} = regexprep(frac_temp{02},' ' ,'');
                CenterCoord_temp                = strsplit(frac_temp{02}, { '[' , ';' , ']' });
                CenterCoord_temp{2}             = strsplit(CenterCoord_temp{2}, { '*' , '/' });
                CenterCoord_temp{3}             = strsplit(CenterCoord_temp{3}, { '*' , '/' });
                CenterCoord_temp{4}             = strsplit(CenterCoord_temp{4}, { '*' , '/' });
                if sum( strcmp('LX' , CenterCoord_temp{2}) ) > 0
                    LX_Index = find( strcmp('LX' , CenterCoord_temp{2}) == 1);
                    CenterCoord_temp{2}{LX_Index} = [];
                    CenterCoord_temp{2}{LX_Index} = num2str(inputMatrix{1}{obj.size+1});
                end
                if length(CenterCoord_temp{2}) == 1
                    obj.FracInput(f).CenterCoord(1,1) = str2double(CenterCoord_temp{2});
                else
                    obj.FracInput(f).CenterCoord(1,1) = str2double(CenterCoord_temp{2}{1}) * str2double(CenterCoord_temp{2}{2}) / str2double(CenterCoord_temp{2}{3});
                end
                if sum( strcmp('LY' , CenterCoord_temp{3}) ) > 0
                    LY_Index = find( strcmp('LY' , CenterCoord_temp{3}) == 1);
                    CenterCoord_temp{3}{LY_Index} = [];
                    CenterCoord_temp{3}{LY_Index} = num2str(inputMatrix{1}{obj.size+2});
                end
                if length(CenterCoord_temp{3}) == 1
                    obj.FracInput(f).CenterCoord(2,1) = str2double(CenterCoord_temp{3});
                else
                    obj.FracInput(f).CenterCoord(2,1) = str2double(CenterCoord_temp{3}{1}) * str2double(CenterCoord_temp{3}{2}) / str2double(CenterCoord_temp{3}{3});
                end
                if sum( strcmp('LZ' , CenterCoord_temp{4}) ) > 0
                    LZ_Index = find( strcmp('LZ' , CenterCoord_temp{4}) == 1);
                    CenterCoord_temp{4}{LZ_Index} = [];
                    CenterCoord_temp{4}{LZ_Index} = num2str(inputMatrix{1}{obj.size+3});
                end
                if length(CenterCoord_temp{4}) == 1
                    obj.FracInput(f).CenterCoord(3,1) = str2double(CenterCoord_temp{4});
                else
                    obj.FracInput(f).CenterCoord(3,1) = str2double(CenterCoord_temp{4}{1}) * str2double(CenterCoord_temp{4}{2}) / str2double(CenterCoord_temp{4}{3});
                end
                %%% Length of the fracture plate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Length_temp = strsplit(frac_temp{03}, { '*' , '/' });
                if sum( strcmp('LX' , Length_temp) ) > 0
                    LX_Index = find( strcmp('LX' , Length_temp) == 1);
                    Length_temp{LX_Index} = [];
                    Length_temp{LX_Index} = num2str(inputMatrix{1}{obj.size+1});
                end
                if length(Length_temp) == 1
                    obj.FracInput(f).Length = str2double(Length_temp);
                else
                    obj.FracInput(f).Length = str2double(Length_temp{1}) * str2double(Length_temp{2}) / str2double(Length_temp{3});
                end
                if sum( strcmp('LY' , Length_temp) ) > 0
                    LY_Index = find( strcmp('LY' , Length_temp) == 1);
                    Length_temp{LY_Index} = [];
                    Length_temp{LY_Index} = num2str(inputMatrix{1}{obj.size+2});
                end
                if length(Length_temp) == 1
                    obj.FracInput(f).Length = str2double(Length_temp);
                else
                    obj.FracInput(f).Length = str2double(Length_temp{1}) * str2double(Length_temp{2}) / str2double(Length_temp{3});
                end
                if sum( strcmp('LZ' , Length_temp) ) > 0
                    LZ_Index = find( strcmp('LZ' , Length_temp) == 1);
                    Length_temp{LZ_Index} = [];
                    Length_temp{LZ_Index} = num2str(inputMatrix{1}{obj.size+3});
                end
                if length(Length_temp) == 1
                    obj.FracInput(f).Length = str2double(Length_temp);
                else
                    obj.FracInput(f).Length = str2double(Length_temp{1}) * str2double(Length_temp{2}) / str2double(Length_temp{3});
                end
                %%% Width of the fracture plate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Width_temp = strsplit(frac_temp{04}, { '*' , '/' });
                if sum( strcmp('LX' , Width_temp) ) > 0
                    LX_Index = find( strcmp('LX' , Width_temp) == 1);
                    Width_temp{LX_Index} = [];
                    Width_temp{LX_Index} = num2str(inputMatrix{1}{obj.size+1});
                end
                if length(Width_temp) == 1
                    obj.FracInput(f).Width = str2double(Width_temp);
                else
                    obj.FracInput(f).Width = str2double(Width_temp{1}) * str2double(Width_temp{2}) / str2double(Width_temp{3});
                end
                if sum( strcmp('LY' , Width_temp) ) > 0
                    LY_Index = find( strcmp('LY' , Width_temp) == 1);
                    Width_temp{LY_Index} = [];
                    Width_temp{LY_Index} = num2str(inputMatrix{1}{obj.size+2});
                end
                if length(Width_temp) == 1
                    obj.FracInput(f).Width = str2double(Width_temp);
                else
                    obj.FracInput(f).Width = str2double(Width_temp{1}) * str2double(Width_temp{2}) / str2double(Width_temp{3});
                end
                if sum( strcmp('LZ' , Width_temp) ) > 0
                    LZ_Index = find( strcmp('LZ' , Width_temp) == 1);
                    Width_temp{LZ_Index} = [];
                    Width_temp{LZ_Index} = num2str(inputMatrix{1}{obj.size+3});
                end
                if length(Width_temp) == 1
                    obj.FracInput(f).Width = str2double(Width_temp);
                else
                    obj.FracInput(f).Width = str2double(Width_temp{1}) * str2double(Width_temp{2}) / str2double(Width_temp{3});
                end
                %%% Rotation angle along Z %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                obj.FracInput(f).RotAngleAlongZ = str2double( frac_temp{05} );
                %%% Rotation angle along central width %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                obj.FracInput(f).RotAngleAlongW = str2double( frac_temp{06} );
                %%% Rotation angle along central length %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                obj.FracInput(f).RotAngleAlongL = str2double( frac_temp{07} );
                %%% Grid resolution ratio between fracture and matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                obj.FracInput(f).GridResRatio   = str2double( frac_temp{08} );
                %%% Grid numbers along the length of the fracture plate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                obj.FracInput(f).GridNumAlongL  = str2double( frac_temp{09} );
                %%% Grid numbers along the width of the fracture plate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                obj.FracInput(f).GridNumAlongW  = str2double( frac_temp{10} );
                %%% The aperture of the fracture plate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                obj.FracInput(f).Aperture       = str2double( frac_temp{11} );
                %%% The porosity of the fracture plate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                obj.FracInput(f).Porosity       = str2double( frac_temp{12} );
                %%% The permeability of the fracture plate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                obj.FracInput(f).Permeability   = str2double( frac_temp{13} );
                %%% ADM Configuration of the fracture %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                frac_temp{14} = regexprep(frac_temp{14},' ' ,'');
                ADM_temp = strsplit(frac_temp{14}, { '[' , ',' , ']' });
                obj.FracInput(f).ADM = [ str2double(ADM_temp{2}) , str2double(ADM_temp{3}) , str2double(ADM_temp{4}) , str2double(ADM_temp{5}) ];
                %%% The coordinates of 3 corners of fracture plate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                frac_temp{15}                 = regexprep(frac_temp{15},' ' ,'');
                CornerCoord_temp              = strsplit(frac_temp{15}, { '[' , ';' , ',' , ']' });
                obj.FracInput(f).CornerCoords = [ str2double(CornerCoord_temp{02}) , str2double(CornerCoord_temp{03}) , str2double(CornerCoord_temp{04})
                                                  str2double(CornerCoord_temp{05}) , str2double(CornerCoord_temp{06}) , str2double(CornerCoord_temp{07})
                                                  str2double(CornerCoord_temp{08}) , str2double(CornerCoord_temp{09}) , str2double(CornerCoord_temp{10}) ];
            end
        end
        
        function Simulation = BuildSimulation(obj, inputMatrix)
            Simulation = simulation_mousa();
            Simulation.Reservoir = obj.BuildReservoir(inputMatrix);
            Simulation.Fractures = obj.BuildFractures(inputMatrix);
            Simulation.Domain = inputMatrix{1}{obj.dim+1};
        end
        
        function Reservoir = BuildReservoir(obj, inputMatrix)
            Reservoir = reservoir_mousa();
            Reservoir.LX = str2double( inputMatrix{1}{obj.size+1} );
            Reservoir.LY = str2double( inputMatrix{1}{obj.size+2} );
            Reservoir.LZ = str2double( inputMatrix{1}{obj.size+3} );
            Reservoir.NX = str2double( inputMatrix{1}{obj.grid+1} );
            Reservoir.NY = str2double( inputMatrix{1}{obj.grid+2} );
            Reservoir.NZ = str2double( inputMatrix{1}{obj.grid+3} );
            CompareAccuracy = str2double( inputMatrix{1}{obj.compareAccuracy+1} );
            Dimension_Type  = inputMatrix{1}{obj.dim+1};
            Reservoir.Discretize(CompareAccuracy,Dimension_Type);
        end
        
        function Fractures = BuildFractures(obj, inputMatrixLX)
            Fractures = fractures_mousa();       
        end
    end
end