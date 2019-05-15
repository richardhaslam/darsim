%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FracGen_Writer  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Mousa HosseiniMehr, MSc Petroleum Engineering, CEG Faculty, TU Delft
% Project: 3D EDFM Package for F-ADM, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build Date : 2017-03-17
% Modified on: 2019-02-11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef writer_FracGen < handle
    properties
        Directory
        File
        VTK_Plotter
    end
    methods
        function obj = writer_FracGen(directory,file)
            %%
            obj.Directory = directory;
            obj.File = file;
        end
        function PrintInfo(obj)
            %%
            disp( '******************* Writing the data into output text file *********************' );
            disp(['Deleting file ', obj.File]);
            delete( strcat(obj.Directory,'/',obj.File) );
            disp(['Writing into file ', obj.File]);
        end
        function WriteOutputTextFile(obj,Geometry)
            ReservoirGrid = Geometry.ReservoirGrid;
            Fracture = Geometry.FracturesGrid.Fracture;
            %% Header
            fid = fopen( strcat(obj.Directory,'/',obj.File) , 'a+' );
            
            fprintf(fid, '%%%% Input file of Fracture properties with connectivities for %s-%s Simulator:\n', Geometry.Domain, Geometry.Type );
            fprintf(fid, '%%%% Convention:\n'                                                        );
            fprintf(fid, '%%%% A %% marks the beginning of a comment.\n'                             );
            fprintf(fid, '%%%% Once started, comments last for the rest of their containing line.\n' );
            fprintf(fid, '%%%% A data block ends with "/".\n'                                        );
            fprintf(fid, '%%%% Keywords are in capital letters.\n'                                   );
            fprintf(fid, '%%%% The index counters follow C++ routine and start from "0". If you use Matlab, please add "1" to every index.\n');
            fprintf(fid, '\n' );
            
            fprintf(fid, '%% The dimension of reservoir:\n');
            fprintf(fid, 'DIMENSION        %8.3f x %8.3f x %8.3f [m^3]\n\n' , ReservoirGrid.LX , ReservoirGrid.LY , ReservoirGrid.LZ );
            
            fprintf(fid, '%% The grid resolution of reservoir:\n');
            fprintf(fid, 'RESERVOIR_GRID   %8.0f x %8.0f x %8.0f [ - ]\n\n' , ReservoirGrid.NX , ReservoirGrid.NY , ReservoirGrid.NZ );
            
            fprintf(fid, '%% Number of fractures:\n');
            fprintf(fid, 'NUM_FRACS     %2.0f\n\n'                     , length(Fracture) );
            
            fprintf(fid, '%% Type of connectivities (either EDFM or pEDFM):\n');
            fprintf(fid, 'TYPE     %s\n', Geometry.Type );
            fprintf(fid, '/\n\n\n' );
            
            %% %% Looping over fractures
            for f = 1 : length(Fracture)
                fprintf(fid, 'FRACTURE\n' );
                fprintf(fid, '%%           fracID   Length      Width       Aperture     Porosity    Permeability   Gridding    ADM_Config [Activation(0,1),Level,Coarsening_ratio_along_length,Coarsening_ratio_along_width]      CornerPoints [PointA] , [PointB] , [PointC] , [PointD]\n');
                fprintf(fid, 'PROPERTIES %2.0f        %1.3e   %1.3e   %1.3e    %1.3e   %1.4e     %1.0fx%1.0f        [%1.0f,%1.0f,%1.0f,%1.0f]                                                                                          [%1.5e;%1.5e;%1.5e] , [%1.5e;%1.5e;%1.5e] , [%1.5e;%1.5e;%1.5e] , [%1.5e;%1.5e;%1.5e]\n' , ...
                    f-1, Fracture(f).Length_AB , Fracture(f).Width_AD  , Fracture(f).Aperture   , Fracture(f).Porosity , ...
                    Fracture(f).Permeability , Fracture(f).N_Length_AB , Fracture(f).N_Width_AD , ...
                    Fracture(f).ADM(1)       , Fracture(f).ADM(2)      , Fracture(f).ADM(3)     , Fracture(f).ADM(4) , ...
                    Fracture(f).PointA(1)    , Fracture(f).PointA(2)   , Fracture(f).PointA(3)  , ...
                    Fracture(f).PointB(1)    , Fracture(f).PointB(2)   , Fracture(f).PointB(3)  , ...
                    Fracture(f).PointC(1)    , Fracture(f).PointC(2)   , Fracture(f).PointC(3)  , ...
                    Fracture(f).PointD(1)    , Fracture(f).PointD(2)   , Fracture(f).PointD(3)  );
                fprintf(fid, '\n' );
                
                %% The Coordinates of grid points
                fprintf(fid, '%% Coordinates of grid points arranged row-wise (through length) then column-wise (through columns)\n');
                fprintf(fid, 'GRID_COORDS_X');
                for j_f = 1 : Fracture(f).N_Width_AD+1
                    fprintf(fid, ' %1.5e', Fracture(f).GridCoords(:,j_f,1)' );
                end
                fprintf(fid, '\n' );
                fprintf(fid, 'GRID_COORDS_Y');
                for j_f = 1 : Fracture(f).N_Width_AD+1
                    fprintf(fid, ' %1.5e', Fracture(f).GridCoords(:,j_f,2)' );
                end
                fprintf(fid, '\n' );
                fprintf(fid, 'GRID_COORDS_Z');
                for j_f = 1 : Fracture(f).N_Width_AD+1
                    fprintf(fid, ' %1.5e', Fracture(f).GridCoords(:,j_f,3)' );
                end
                fprintf(fid, '\n\n' );
                
                %% Looping over each fracture's grid cells
                for If = 1 : Fracture(f).N_Total
                    if ~strcmp(Geometry.Type,'pEDFM')
                        Fracture(f).CI_rock_pEDFM{If}=[];
                        Fracture(f).NumOf_fracConn_pEDFM(If)=0;
                    end
                    fprintf(fid, '%%                    fracCellID         #rock EDFM    #rock pEDFM    #frac EDFM    #frac pEDFM\n' );
                    fprintf(fid, 'FRACCELL                %3.0f             %3.0f            %3.0f           %3.0f           %3.0f\n' , ...
                        If-1 , size(Fracture(f).CI_rock_EDFM{If},1) , size(Fracture(f).CI_rock_pEDFM{If},1) , Fracture(f).NumOf_fracConn_EDFM(If), Fracture(f).NumOf_fracConn_pEDFM(If) );
                    
                    % Writing Rock EDFM Connectivities
                    fprintf(fid, '%%                    rockCellID           CI_rock_EDFM\n' );
                    for im = 1 : size( Fracture(f).CI_rock_EDFM{If} , 1 )
                        fprintf(fid, 'ROCK_CONN_EDFM        %5.0f               %1.5e\n' , Fracture(f).CI_rock_EDFM{If}(im,1)-1 , Fracture(f).CI_rock_EDFM{If}(im,2) );
                    end
                    
                    % Writing Rock pEDFM Connectivities
                    if strcmp(Geometry.Type,'pEDFM')
                        fprintf(fid, '%%                    rockCellID           CI_rock_pEDFM\n' );
                        for im = 1 : size( Fracture(f).CI_rock_pEDFM{If} , 1 )
                            fprintf(fid, 'ROCK_CONN_pEDFM       %5.0f               %1.5e\n' , Fracture(f).CI_rock_pEDFM{If}(im,1)-1 , Fracture(f).CI_rock_pEDFM{If}(im,2) );
                        end
                    end
                    
                    % Writing Frac EDFM Connectivities
                    if Fracture(f).NumOf_fracConn_EDFM(If) > 0
                        fprintf(fid, '%%                fracID and fracCellID    CI_frac_EDFM\n' );
                        for g = 1 : length(Fracture)
                            if ( ~isempty(Fracture(f).CI_frac_EDFM{g}) )
                                if ( ~isempty(Fracture(f).CI_frac_EDFM{g}{If}) )
                                    for ig = 1 : size ( Fracture(f).CI_frac_EDFM{g}{If} , 1 )
                                        if ( ~isempty(Fracture(f).CI_frac_EDFM{g}{If}(ig,1)) )
                                            fprintf(fid, 'FRAC_CONN_EDFM   %5.0f    %5.0f           %1.5e\n' , g-1 , Fracture(f).CI_frac_EDFM{g}{If}(ig,1)-1 , Fracture(f).CI_frac_EDFM{g}{If}(ig,2) );
                                        end
                                    end
                                end
                            end
                        end
                    end
                    
                    % Writing Frac pEDFM Connectivities
                    if strcmp(Geometry.Type,'pEDFM')
                        if Fracture(f).NumOf_fracConn_pEDFM(If) > 0
                            fprintf(fid, '%%                fracID and fracCellID    CI_frac_pEDFM\n' );
                            for g = 1 : length(Fracture)
                                if ( ~isempty(Fracture(f).CI_frac_pEDFM{g}) )
                                    if ( ~isempty(Fracture(f).CI_frac_pEDFM{g}{If}) )
                                        for ig = 1 : size ( Fracture(f).CI_frac_pEDFM{g}{If} , 1 )
                                            if ( ~isempty(Fracture(f).CI_frac_pEDFM{g}{If}(ig,1)) )
                                                fprintf(fid, 'FRAC_CONN_pEDFM  %5.0f    %5.0f           %1.5e\n' , g-1 , Fracture(f).CI_frac_pEDFM{g}{If}(ig,1)-1 , Fracture(f).CI_frac_pEDFM{g}{If}(ig,2) );
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                    
                    fprintf(fid, '\n' );
                end
                fprintf(fid, '/\n\n\n' );
            end
            
            %% The pEDFM alpha correction for reservoir tranmissibilities:
            if strcmp(Geometry.Type,'pEDFM')
                fprintf(fid, '%%%% The pEDFM alpha correction for reservoir tranmissibilities:\n' );
                fprintf(fid, '%%%% The transmissiblities in the simulation software should be multiplied by "1-alpha".\n\n' );
                fprintf(fid, '%% Tx( NX+1,NY,NZ )        i       j       k       alpha\n' );
                for k = 1:ReservoirGrid.NZ
                    for j = 1:ReservoirGrid.NY
                        for i = 1:ReservoirGrid.NX+1
                            if ReservoirGrid.alpha_Tx(i,j,k) > ReservoirGrid.Epsilon
                                fprintf(fid, 'ROCK_ALPHA_TX           %3.0f     %3.0f     %3.0f       %1.5e\n' , i-1,j-1,k-1  , ReservoirGrid.alpha_Tx(i,j,k) );
                            end
                        end
                    end
                end
                fprintf(fid, '\n' );
                
                fprintf(fid, '%% Ty( NX,NY+1,NZ )        i       j       k       alpha\n' );
                for k = 1:ReservoirGrid.NZ
                    for j = 1:ReservoirGrid.NY+1
                        for i = 1:ReservoirGrid.NX
                            if ReservoirGrid.alpha_Ty(i,j,k) > ReservoirGrid.Epsilon
                                fprintf(fid, 'ROCK_ALPHA_TY           %3.0f     %3.0f     %3.0f       %1.5e\n' , i-1,j-1,k-1  , ReservoirGrid.alpha_Ty(i,j,k) );
                            end
                        end
                    end
                end
                fprintf(fid, '\n' );
                
                fprintf(fid, '%% Tz( NX,NY,NZ+1 )        i       j       k       alpha\n' );
                for k = 1:ReservoirGrid.NZ+1
                    for j = 1:ReservoirGrid.NY
                        for i = 1:ReservoirGrid.NX
                            if ReservoirGrid.alpha_Tz(i,j,k) > ReservoirGrid.Epsilon
                                fprintf(fid, 'ROCK_ALPHA_TZ           %3.0f     %3.0f     %3.0f       %1.5e\n' , i-1,j-1,k-1  , ReservoirGrid.alpha_Tz(i,j,k) );
                            end
                        end
                    end
                end
                fprintf(fid, '\n\n\n' );
            end
            
            %% The pEDFM alpha correction for fractures tranmissibilities:
            if strcmp(Geometry.Type,'pEDFM')
                fprintf(fid, '%%%% The pEDFM alpha correction for fractures tranmissibilities:\n' );
                fprintf(fid, '%%%% The transmissiblities in the simulation software should be multiplied by "1-alpha".\n\n' );
                fprintf(fid, '%% Tx( NX+1,NY )           frac#   i_f     j_f     alpha\n' );
                for f=1:length(Fracture)
                    for j_f = 1 : Fracture(f).N_Width_AD
                        for i_f = 1 : Fracture(f).N_Length_AB + 1
                            if Fracture(f).alpha_Tx(i_f,j_f) > ReservoirGrid.Epsilon
                                fprintf(fid, 'FRAC_ALPHA_TX           %3.0f     %3.0f     %3.0f       %1.5e\n' , f-1,i_f-1,j_f-1  , Fracture(f).alpha_Tx(i_f,j_f) );
                            end
                        end
                    end
                end
                fprintf(fid, '\n' );
                
                fprintf(fid, '%% Ty( NX,NY+1 )           frac#   i_f     j_f     alpha\n' );
                for f=1:length(Fracture)
                    for j_f = 1 : Fracture(f).N_Width_AD + 1
                        for i_f = 1 : Fracture(f).N_Length_AB
                            if Fracture(f).alpha_Ty(i_f,j_f) > ReservoirGrid.Epsilon
                                fprintf(fid, 'FRAC_ALPHA_TY           %3.0f     %3.0f     %3.0f       %1.5e\n' , f-1,i_f-1,j_f-1  , Fracture(f).alpha_Ty(i_f,j_f) );
                            end
                        end
                    end
                end
            end
            
            %%
            fclose(fid);
            disp('The output file is ready!');
            
        end
    end
end