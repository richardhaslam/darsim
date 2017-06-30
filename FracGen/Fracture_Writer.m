%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Fracture_Writer  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Mousa HosseiniMehr, MSc Petroleum Engineering, CEG Faculty, TU Delft
% Project: 3D EDFM Package for F-ADM, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build Date : 2017-03-17
% Modified on: 2017-03-17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Fracture_Writer(File, obj)
disp( '******************* Writing the data into output text file *********************' );
disp(['Deleting file ', File]);
delete( File );
disp(['Writing into file ', File]);

%%
fid = fopen( File , 'a+' );

fprintf(fid, ['%%%% Input file of Fracture properties with connectivities for ', obj.Simulation.Domain, '-pEDFM Simulator\n'] );
fprintf(fid, '%%%% Convention:\n'                                                            );
fprintf(fid, '%%%% A %% marks the beginning of a comment\n'                                  );
fprintf(fid, '%%%% Once started, comments last for the rest of their containing line\n'      );
fprintf(fid, '%%%% A data block ends with /\n'                                               );
fprintf(fid, '%%%% Keywords are in capital letters\n'                                        );
fprintf(fid, '\n' );

if obj.Simulation.Domain == '2D'
    %% 2D Case
    fprintf(fid, '%%Dimension   %8.1f x %8.1f [m^2]\n' , obj.Simulation.Reservoir.LX , obj.Simulation.Reservoir.LY );
    fprintf(fid, '%%Grid        %8.0f x %8.0f [ - ]\n' , obj.Simulation.Reservoir.NX , obj.Simulation.Reservoir.NY );
    fprintf(fid, 'NUM_FRACS     %2.0f\n' , length(obj.Simulation.Fractures)             );
    fprintf(fid, '\n' );
    fprintf(fid, '/\n\n\n' );
    for f = 1 : length(obj.Simulation.Fractures)
        fprintf(fid, 'FRACTURE\n' );
        fprintf(fid, '%%        	fracID      Length    Aperture      Porosity      Permeability      Gridding      CornerPoints [PointA] , [PointB]\n');
        fprintf(fid, 'PROPERTIES    %2.0f       %1.5e     %1.5e         %1.5e         %1.5e             %1.0f         [%1.5e;%1.5e] , [%1.5e;%1.5e]\n' , ...
                f-1 , obj.Simulation.Fractures(f).Length_AB , obj.Simulation.Fractures(f).Aperture , obj.Simulation.Fractures(f).Porosity , ...
                obj.Simulation.Fractures(f).Permeability , obj.Simulation.Fractures(f).N_Length_AB , ...
                     obj.Simulation.Fractures(f).PointA(1) , obj.Simulation.Fractures(f).PointA(2) , ...
                     obj.Simulation.Fractures(f).PointB(1) , obj.Simulation.Fractures(f).PointB(2) );
                 
        for i_f = 1 : obj.Simulation.Fractures(f).N_Length_AB
            
            fprintf(fid, '%%        fracCellID	#perf	#proj\n' );
            fprintf(fid, 'FRACCELL	%2.0f	%2.0f	%2.0f	%2.0f\n' , i_f-1 , size(obj.Simulation.Fractures(f).areaFrac_matCell{i_f},1) , 0 );

            fprintf(fid, '%%            matCellID	area	aveDist\n' );
            for im = 1 : size( obj.Simulation.Fractures(f).areaFrac_matCell{i_f} , 1 )
                fprintf(fid, 'ROCK_CONN	%2.0f	%1.5e	%1.5e\n' , obj.Simulation.Fractures(f).areaFrac_matCell{i_f}{im,1}-1 , obj.Simulation.Fractures(f).areaFrac_matCell{i_f}{im,2} , obj.Simulation.Fractures(f).aveDist_matCell{i_f}{im,2} );
            end
      
        end
        
        fprintf(fid, '%%           [neighboring cellIDs in this fracture] , neighboring fracID , [neighboring cellIDs in the other fracture]\n' );
        for g = 1 : length(obj.Simulation.Fractures)
            if ~isempty( obj.Simulation.Fractures(f).StarDelta2D_Neighbor{g,1} )
                
                fprintf(fid, 'FRAC_CONN	[%1.0f,%1.0f] , %1.0f , [%1.0f,%1.0f]\n' , obj.Simulation.Fractures(f).StarDelta2D_Neighbor{g,1}.itself(1) , obj.Simulation.Fractures(f).StarDelta2D_Neighbor{g,1}.itself(2) , ...
                                                                             g-1 , obj.Simulation.Fractures(f).StarDelta2D_Neighbor{g,1}.other(1)  , obj.Simulation.Fractures(f).StarDelta2D_Neighbor{g,1}.other(2)  );
            end
        end

    end


elseif obj.Simulation.Domain == '3D'
    %% 3D CaSE    
    fprintf(fid, '%%Dimension   %8.1f x %8.1f x %8.1f [m^3]\n' , obj.Simulation.Reservoir.LX , obj.Simulation.Reservoir.LY , obj.Simulation.Reservoir.LZ );
    fprintf(fid, '%%Grid        %8.0f x %8.0f x %8.0f [ - ]\n' , obj.Simulation.Reservoir.NX , obj.Simulation.Reservoir.NY , obj.Simulation.Reservoir.NZ );
    fprintf(fid, 'NUM_FRACS     %2.0f\n' , length(obj.Simulation.Fractures)             );
    fprintf(fid, '\n' );
    fprintf(fid, '/\n\n\n' );
    for f = 1 : length(obj.Simulation.Fractures)
        fprintf(fid, 'FRACTURE\n' );
        fprintf(fid, '%%        	fracID      Length      Width      Aperture      Porosity      Permeability      Gridding      CornerPoints [PointA] , [PointB] , [PointC] , [PointD]\n');
        fprintf(fid, 'PROPERTIES    %2.0f       %1.5e       %1.5e      %1.5e         %1.5e         %1.5e             %1.0fx%1.0f   [%1.5e;%1.5e;%1.5e] , [%1.5e;%1.5e;%1.5e] , [%1.5e;%1.5e;%1.5e] , [%1.5e;%1.5e;%1.5e]\n' , ...
                f-1 , obj.Simulation.Fractures(f).Length_AB , obj.Simulation.Fractures(f).Width_AD , obj.Simulation.Fractures(f).Aperture , obj.Simulation.Fractures(f).Porosity , ...
                obj.Simulation.Fractures(f).Permeability , obj.Simulation.Fractures(f).N_Length_AB , obj.Simulation.Fractures(f).N_Width_AD , ...
                     obj.Simulation.Fractures(f).PointA(1) , obj.Simulation.Fractures(f).PointA(2) , obj.Simulation.Fractures(f).PointA(3) , ...
                     obj.Simulation.Fractures(f).PointB(1) , obj.Simulation.Fractures(f).PointB(2) , obj.Simulation.Fractures(f).PointB(3) , ...
                     obj.Simulation.Fractures(f).PointC(1) , obj.Simulation.Fractures(f).PointC(2) , obj.Simulation.Fractures(f).PointC(3) , ...
                     obj.Simulation.Fractures(f).PointD(1) , obj.Simulation.Fractures(f).PointD(2) , obj.Simulation.Fractures(f).PointD(3) );
                 
        fprintf(fid, '%% Coordinates of grid points arranged row-wise (through length) then column-wise (through columns)\n');
        fprintf(fid, 'GRID_COORDS_X');
        for j_f = 1 : obj.Simulation.Fractures(f).N_Width_AD+1
            fprintf(fid, ' %1.5e', obj.Simulation.Fractures(f).GridCoords(:,j_f,1)' );
        end
        fprintf(fid, '\n' );
        fprintf(fid, 'GRID_COORDS_Y');
        for j_f = 1 : obj.Simulation.Fractures(f).N_Width_AD+1
            fprintf(fid, ' %1.5e', obj.Simulation.Fractures(f).GridCoords(:,j_f,2)' );
        end
        fprintf(fid, '\n' );
        fprintf(fid, 'GRID_COORDS_Z');
        for j_f = 1 : obj.Simulation.Fractures(f).N_Width_AD+1
            fprintf(fid, ' %1.5e', obj.Simulation.Fractures(f).GridCoords(:,j_f,3)' );
        end
        fprintf(fid, '\n' );
        

        for j_f = 1 : obj.Simulation.Fractures(f).N_Width_AD
            for i_f = 1 : obj.Simulation.Fractures(f).N_Length_AB
                If = Index_Fracture_2D( obj.Simulation.Fractures(f).N_Length_AB , obj.Simulation.Fractures(f).N_Width_AD , i_f , j_f );
                fprintf(fid, '%%        fracCellID	#perf	#proj	#fracConn\n' );
                fprintf(fid, 'FRACCELL	%2.0f	%2.0f	%2.0f	%2.0f\n' , If-1 , size(obj.Simulation.Fractures(f).areaFrac_matCell{If},1) , 0 , obj.Simulation.Fractures(f).NumOf_fracCellConn(If) );

                fprintf(fid, '%%            matCellID	area	aveDist\n' );
                for im = 1 : size( obj.Simulation.Fractures(f).areaFrac_matCell{If} , 1 )
                    fprintf(fid, 'ROCK_CONN	%2.0f	%1.5e	%1.5e\n' , obj.Simulation.Fractures(f).areaFrac_matCell{If}{im,1}-1 , obj.Simulation.Fractures(f).areaFrac_matCell{If}{im,2} , obj.Simulation.Fractures(f).aveDist_matCell{If}{im,2} );
                end

                fprintf(fid, '%%            fracID	cellID	area	aveDist\n' );
                for g = 1 : length(obj.Simulation.Fractures)
                   if ( ~isempty( obj.Simulation.Fractures(f).areaFrac_fracCell{g} ) )
                       if ( ~isempty( obj.Simulation.Fractures(f).areaFrac_fracCell{g}{If} ) )
                           for ig = 1 : size ( obj.Simulation.Fractures(f).areaFrac_fracCell{g}{If} , 1 )
                               if ( ~isempty( obj.Simulation.Fractures(f).areaFrac_fracCell{g}{If}{ig,1} ) )
                                   fprintf(fid, 'FRAC_CONN	%2.0f	%3.0f	%1.5e	%1.5e\n' , g-1 , obj.Simulation.Fractures(f).areaFrac_fracCell{g}{If}{ig,1}-1 , obj.Simulation.Fractures(f).areaFrac_fracCell{g}{If}{ig,2} , obj.Simulation.Fractures(f).aveDist_fracCell{g}{If}{ig,2} );
                               end
                           end
                       end
                   end
                end

            end
        end
        fprintf(fid, '/\n\n\n' );
    end
    
else
    
    Error( 'The Domain (2D or 3D) is not mentioned correctly in the input file! Is there a typo?' );
    
end

fclose(fid);

disp('The output file is ready!');

end