% Simulation class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Mousa HosseiniMehr
%TU Delft
%Created: 12 July 2016
%Last modified: 17 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef geometry_FracGen
    properties
        Type
        Domain
        ReservoirGrid
        FracturesGrid
    end
    methods
        function PrintInfo(obj)
            disp( ['************* Generating ', obj.Type, ' Fractures for ', obj.Domain,' Domain Simulation **************']);
        end
        function Compute_EDFM_Connectivities(obj)
            %%
            Fracture = obj.FracturesGrid.Fracture;
            Epsilon = obj.ReservoirGrid.Epsilon;

            %% Intializing necessary variables
            for f = 1 : length(Fracture)
                % Coordinates of intersections between each two fracture plates:
                Fracture(f).intersectCoord_fracturePlate = cell( length(Fracture) , 1 );
                % Overlap status of fracture f cells due to intersection with fracture g:
                Fracture(f).Overlap_frac = cell( length(Fracture) , 1 );
                
                Fracture(f).intersectCoord_rock = cell( Fracture(f).N_Total , 1 );   % Coordinates of intersections between each fracture cell and each matrix cell
                Fracture(f).      Af_rock_EDFM  = cell( Fracture(f).N_Total , 1 );   % Area fraction of each fracture cell inside each matrix cube
                Fracture(f). avgDist_rock_EDFM  = cell( Fracture(f).N_Total , 1 );   % Average distance between each fracture cell and each matrix cube
                Fracture(f).      CI_rock_EDFM  = cell( Fracture(f).N_Total , 1 );   % Geometrical Transmissibility between each fracture cell and each matrix cube
                Fracture(f).  Af_rock_EDFM_sum  = 0;                                 % Summation of all area fractions (only to validate)
                Fracture(f).intersectCoord_frac = cell( length(Fracture) , 1 );      % Coordinates of intersections between each two fracture cells of distinct fracture plates
                Fracture(f).       Af_frac_EDFM = cell( length(Fracture) , 1 );      % Area fraction of of intersection line between each two fracture cells of distinct fracture plates
                Fracture(f).  avgDist_frac_EDFM = cell( length(Fracture) , 1 );      % Average distance between each fracture cell and the intersection line between each two fracture cells of distinct fracture plates
                Fracture(f).       CI_frac_EDFM = cell( length(Fracture) , 1 );      % Geometrical Transmissibility between each two non-neighboring fracture cells of distinct fracture plates
                Fracture(f).NumOf_fracConn_EDFM = zeros( Fracture(f).N_Total , 1);   % Number of non-neighboring fracture cells connections
                Fracture(f).      Af_rock_pEDFM = cell(Fracture(f).N_Total,1);
                Fracture(f). avgDist_rock_pEDFM = cell(Fracture(f).N_Total,1);
                Fracture(f).      CI_rock_pEDFM = cell(Fracture(f).N_Total,1);
                Fracture(f).      Af_frac_pEDFM = cell(length(Fracture),1);
                Fracture(f). avgDist_frac_pEDFM = cell(length(Fracture),1);
                Fracture(f).      CI_frac_pEDFM = cell(length(Fracture),1);
                Fracture(f).NumOf_fracConn_pEDFM = zeros( Fracture(f).N_Total , 1);
            end
            
            %% Intersections of fracture plates with each other (if any)
            fprintf('Obtaining the intersection of fracture plates ...\n');
            for f = 1 : length(Fracture)
                for g = f+1 : length(Fracture)
                    
                    % Obtaining the intersection Points
                    [Geostatus, intersectPoints] = Fracture(f).Obtain_PlaneSegment_PlaneSegment_Intersection( Fracture(g), Epsilon );
                    intersectCoord_fracturePlate_final = [];
                    
                    if ~isempty(intersectPoints) && (Geostatus.areCoplanar ~= 1 )
                        % Processing the intersection Points
                        for nr = 1 : size(intersectPoints,2)
                            
                            % Removing the Points that are not inside the fracture rectangles
                            isInside1 = Fracture(f).Is_Point_Inside_PlaneSegment(intersectPoints(:,nr), Epsilon);
                            isInside2 = Fracture(g).Is_Point_Inside_PlaneSegment(intersectPoints(:,nr), Epsilon);
                            if ( isInside1 == 0 ) || ( isInside2 == 0 )
                                intersectPoints(:,nr) = [ NaN ; NaN ; NaN ];
                                continue;
                            end
                            
                            % Removing the points that are repeated
                            subtracted = ( intersectPoints - intersectPoints(:,nr) );
                            subtracted = sqrt( subtracted(1,:).^2 + subtracted(2,:).^2 + subtracted(3,:).^2 );
                            if ~isempty( find( subtracted([1:nr-1 , nr+1:end]) <= Epsilon , 1 ) )
                                intersectPoints(:,nr) = [ NaN ; NaN ; NaN ];
                                continue;
                            end
                            
                            if ~isnan(intersectPoints(1,nr))
                                intersectCoord_fracturePlate_final = [ intersectCoord_fracturePlate_final , intersectPoints(:,nr) ];
                            end
                        end
                        
                        % If there is no or only one final intersection Point, do not plot
                        if size(intersectCoord_fracturePlate_final,2) < 2,  continue;  end
                        
                        % Writing the final intersection Points into the fracture data structure
                        Fracture(f).intersectCoord_fracturePlate{g} = intersectCoord_fracturePlate_final;
                        Fracture(g).intersectCoord_fracturePlate{f} = intersectCoord_fracturePlate_final;
                        Fracture(f).Overlap_frac{g} = zeros( Fracture(f).N_Length_AB * Fracture(f).N_Width_AD , 1 );
                        Fracture(g).Overlap_frac{f} = zeros( Fracture(g).N_Length_AB * Fracture(g).N_Width_AD , 1 );
                    end
                end
            end
            fprintf('---------------------------------------------------------\n');
            
            %% Obtaining fracture - matrix connectivities
            obj.Compute_EDFM_Fracture_Matrix_Connectivities();
            %% Obtaining fracture - fracture connectivities
            obj.Compute_EDFM_Fracture_Fracture_Connectivities();
        end
        function Compute_EDFM_Fracture_Matrix_Connectivities(obj)
            Fracture = obj.FracturesGrid.Fracture;
            %% Assigning reservoir properties
            % Number of reservoir grid cells in each direction
            NX = obj.ReservoirGrid.NX;
            NY = obj.ReservoirGrid.NY;
            NZ = obj.ReservoirGrid.NZ;
            % Coordinates of reservoir cell centers in each direction
            Xcm = obj.ReservoirGrid.Xcm;
            Ycm = obj.ReservoirGrid.Ycm;
            Zcm = obj.ReservoirGrid.Zcm;
            % Coordinates of reservoir grid nodes in each direction
            Xim = obj.ReservoirGrid.Xim;
            Yim = obj.ReservoirGrid.Yim;
            Zim = obj.ReservoirGrid.Zim;
            %% Intersections of fracture cells with each matrix cell
            fprintf('Obtaining fracture - matrix connectivities:\n');
            fprintf('---> Fracture ');
            for f = 1 : length(Fracture)
                if (f>1),  fprintf(repmat('\b', 1, 5+27));  end
                fprintf('%02d/%02d',f,length(Fracture));
                fprintf(' ---> Grid cell ');
                for j_f = 1 : Fracture(f).N_Width_AD
                    for i_f = 1 : Fracture(f).N_Length_AB
                        If = Index_Local_To_Global_3D( Fracture(f).N_Length_AB , Fracture(f).N_Width_AD , 1 , i_f , j_f , 1 );
                        if (If>1),  fprintf(repmat('\b', 1, 11));  end
                        fprintf('%05d/%05d',If,Fracture(f).N_Length_AB*Fracture(f).N_Width_AD);
                        Index_matIntersect = 0;
                        
                        % The matrix cell closest to fracture cell for the first intersection check
                        [ ~ , i_1st ] = min( ( Fracture(f).CellCenterCoordsV1(i_f,j_f,1) - Xcm ).^2 );
                        [ ~ , j_1st ] = min( ( Fracture(f).CellCenterCoordsV1(i_f,j_f,2) - Ycm ).^2 );
                        [ ~ , k_1st ] = min( ( Fracture(f).CellCenterCoordsV1(i_f,j_f,3) - Zcm ).^2 );
                        Im_List  = Index_Local_To_Global_3D( NX,NY,NZ , i_1st , j_1st , k_1st );
                        Im_count = 1;
                        
                        while Im_count <= length(Im_List)
                            Im = Im_List(Im_count);
                            Im_count = Im_count +1;
                            
                            % Retrieving i,j,k from the "Im" index
                            [i,j,k] = Index_Global_To_Local_3D(NX, NY, NZ, Im);

                            % Corner Points of the fracture cell
                            PointA = [ Fracture(f).GridCoords(i_f  ,j_f  ,1) ; Fracture(f).GridCoords(i_f  ,j_f  ,2) ; Fracture(f).GridCoords(i_f  ,j_f  ,3) ];
                            PointB = [ Fracture(f).GridCoords(i_f  ,j_f+1,1) ; Fracture(f).GridCoords(i_f  ,j_f+1,2) ; Fracture(f).GridCoords(i_f  ,j_f+1,3) ];
                            PointC = [ Fracture(f).GridCoords(i_f+1,j_f+1,1) ; Fracture(f).GridCoords(i_f+1,j_f+1,2) ; Fracture(f).GridCoords(i_f+1,j_f+1,3) ];
                            PointD = [ Fracture(f).GridCoords(i_f+1,j_f  ,1) ; Fracture(f).GridCoords(i_f+1,j_f  ,2) ; Fracture(f).GridCoords(i_f+1,j_f  ,3) ];
                            Plane_fracCell = planeSegment_FracGen(PointA,PointB,PointC,PointD);
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % Obtaining the temporary and finqal intersection points between the matrix cell and fracture cell
                            [intersectCoordFinal, intersectCoordTemp , areCoplanar] = obj.Obtain_PlaneSegment_Cube_Intersection(i,j,k, Plane_fracCell);
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % Add the neighboring matrix cells to the list for intersection check if it is the first try
                            % but no intersection occurs, so another one must be checked
                            if ( ( isempty(intersectCoordTemp) ) || ( size(intersectCoordTemp,2) < 3 ) ) && Im_count == 2
                                Im_next = Index_Local_To_Global_3D( NX,NY,NZ,max(i-1,1 ),j          ,k           );  if ~ismember(Im_next,Im_List),  Im_List = [Im_List;Im_next];  end
                                Im_next = Index_Local_To_Global_3D( NX,NY,NZ,min(i+1,NX),j          ,k           );  if ~ismember(Im_next,Im_List),  Im_List = [Im_List;Im_next];  end
                                Im_next = Index_Local_To_Global_3D( NX,NY,NZ,i          ,max(j-1,1 ),k           );  if ~ismember(Im_next,Im_List),  Im_List = [Im_List;Im_next];  end
                                Im_next = Index_Local_To_Global_3D( NX,NY,NZ,i          ,min(j+1,NY),k           );  if ~ismember(Im_next,Im_List),  Im_List = [Im_List;Im_next];  end
                                Im_next = Index_Local_To_Global_3D( NX,NY,NZ,i          ,j          ,max(k-1,1 ) );  if ~ismember(Im_next,Im_List),  Im_List = [Im_List;Im_next];  end
                                Im_next = Index_Local_To_Global_3D( NX,NY,NZ,i          ,j          ,min(k+1,NZ) );  if ~ismember(Im_next,Im_List),  Im_List = [Im_List;Im_next];  end
                            end
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % Assigning the index of matrix cell to the array
                            if Im_count > 2
                                if ( isempty(intersectCoordTemp) ),  continue;  end
                            end
                            
                            % Add the neighboring cells to the list for intersection check
                            Im_next = Index_Local_To_Global_3D( NX,NY,NZ,max(i-1,1 ),j          ,k           );  if ~ismember(Im_next,Im_List),  Im_List = [Im_List;Im_next];  end
                            Im_next = Index_Local_To_Global_3D( NX,NY,NZ,min(i+1,NX),j          ,k           );  if ~ismember(Im_next,Im_List),  Im_List = [Im_List;Im_next];  end
                            Im_next = Index_Local_To_Global_3D( NX,NY,NZ,i          ,max(j-1,1 ),k           );  if ~ismember(Im_next,Im_List),  Im_List = [Im_List;Im_next];  end
                            Im_next = Index_Local_To_Global_3D( NX,NY,NZ,i          ,min(j+1,NY),k           );  if ~ismember(Im_next,Im_List),  Im_List = [Im_List;Im_next];  end
                            Im_next = Index_Local_To_Global_3D( NX,NY,NZ,i          ,j          ,max(k-1,1 ) );  if ~ismember(Im_next,Im_List),  Im_List = [Im_List;Im_next];  end
                            Im_next = Index_Local_To_Global_3D( NX,NY,NZ,i          ,j          ,min(k+1,NZ) );  if ~ismember(Im_next,Im_List),  Im_List = [Im_List;Im_next];  end

                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % Writing the final intersection Points into the fracture data structure
                            if ( size(intersectCoordFinal,2) < 3 ),  continue;  end
                            
                            Index_matIntersect = Index_matIntersect + 1;
                            Fracture(f).intersectCoord_rock{If}{Index_matIntersect,1} = Im;
                            Fracture(f).       Af_rock_EDFM{If}(Index_matIntersect,1) = Im;
                            Fracture(f).  avgDist_rock_EDFM{If}(Index_matIntersect,1) = Im;
                            Fracture(f).       CI_rock_EDFM{If}(Index_matIntersect,1) = Im;
                            
                            Fracture(f).intersectCoord_rock{If}{Index_matIntersect,2} = intersectCoordFinal;
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % Calculating the area fraction
                            Fracture(f).Af_rock_EDFM{If}(Index_matIntersect,2) = 0;
                            for n = 2 : size( intersectCoordFinal , 2 )  - 1
                                Fracture(f).Af_rock_EDFM{If}(Index_matIntersect,2) = Fracture(f).Af_rock_EDFM{If}(Index_matIntersect,2) + ...
                                    Triangle_Area_3D( intersectCoordFinal(:,1) , intersectCoordFinal(:,n) , intersectCoordFinal(:,n+1) );
                            end
                            
                            % Doubling the area fraction if the cell is not coplanar to any matrix cube edges
                            if areCoplanar == 0
                                Fracture(f).Af_rock_EDFM{If}(Index_matIntersect,2) = Fracture(f).Af_rock_EDFM{If}(Index_matIntersect,2) * 2;
                            end
                            Fracture(f).Af_rock_EDFM_sum = Fracture(f).Af_rock_EDFM_sum + Fracture(f).Af_rock_EDFM{If}(Index_matIntersect,2);
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % Obtaining the average distance between the matrix cube cell and the fracture plate cell
                            Fracture(f).avgDist_rock_EDFM{If}(Index_matIntersect,2) = 0;
                            CubeCornerStart = [ Xim(i)  ;Yim(j)  ;Zim(k)  ];
                            CubeCornerEnd   = [ Xim(i+1);Yim(j+1);Zim(k+1)];
                            Refinement = 6;
                            Fracture(f).avgDist_rock_EDFM{If}(Index_matIntersect,2) = ...
                                Average_Distance_Cube_From_Plane( Fracture(f).Equation , CubeCornerStart , CubeCornerEnd , Refinement );
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % Obtaining the geometrical transmissibility between the matrix cube cell and the fracture plate cell ( CI = Af/<d> )
                            Fracture(f).CI_rock_EDFM{If}(Index_matIntersect,2) = ...
                            Fracture(f).Af_rock_EDFM{If}(Index_matIntersect,2) / ...
                            Fracture(f).avgDist_rock_EDFM{If}(Index_matIntersect,2);
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % Adding the matrix cell index to the list of
                            % all overlapped matrix cells inside the
                            % "reservoir" class and "fracture" class
                            if ~ismember(Im,obj.ReservoirGrid.OvelapList)
                                obj.ReservoirGrid.OvelapList = [obj.ReservoirGrid.OvelapList;Im];
                            end
                            if ~ismember(Im,Fracture(f).MatrixOverlapList)
                                Fracture(f).MatrixOverlapList = [Fracture(f).MatrixOverlapList;Im];
                            end
                            
                        end % End of while-loop
                        
                        % Obtaining the overlapped cells due to the intersection line between fractures f and g
                        max_dist = sqrt( Fracture(f).D_Length_AB^2 + Fracture(f).D_Width_AD^2 )/2;
                        for g = 1 : length(Fracture)
                            if ~isempty(Fracture(f).intersectCoord_fracturePlate{g})
                                p1 = Fracture(f).intersectCoord_fracturePlate{g}(:,1);
                                p2 = Fracture(f).intersectCoord_fracturePlate{g}(:,2);
                                Fracture(f).Overlap_frac{g}((j_f-1)*Fracture(f).N_Length_AB+i_f) = ( norm( cross(Fracture(f).CellCenterCoordsV2((j_f-1)*Fracture(f).N_Length_AB+i_f,:)'-p1 , p2-p1) ) / norm(p2-p1) ) <= max_dist;
                            end
                        end
                        Fracture(f).     Af_rock_EDFM{If} = sortrows(Fracture(f).     Af_rock_EDFM{If});
                        Fracture(f).avgDist_rock_EDFM{If} = sortrows(Fracture(f).avgDist_rock_EDFM{If});
                        Fracture(f).     CI_rock_EDFM{If} = sortrows(Fracture(f).     CI_rock_EDFM{If});
                    end % End of fractre 1st for-loop
                end % End of fractre 2nd for-loop
                Fracture(f).MatrixOverlapList = sort(Fracture(f).MatrixOverlapList);
            end % End of main fracture for-loop
            obj.ReservoirGrid.OvelapList = sort(obj.ReservoirGrid.OvelapList);
            
            fprintf('\n---------------------------------------------------------\n');
        end
        function Compute_EDFM_Fracture_Fracture_Connectivities(obj)
            Fracture = obj.FracturesGrid.Fracture;
            Epsilon = obj.ReservoirGrid.Epsilon;
            %% Intersections of fracture cells of distinct fractures with each other (if any)
            fprintf('Obtaining fracture - fracture connectivities:\n');
            for f = 1 : length(Fracture)
                for g = f+1 : length(Fracture)
                    % Finding the intersections between the cells of distinct fractures
                    if ~isempty(Fracture(f).intersectCoord_fracturePlate{g})
                        fprintf('Fracture %02d <---> %02d\n',f,g);
                        Fracture(f).intersectCoord_frac{g} = cell( Fracture(f).N_Total , 1 );
                        Fracture(g).intersectCoord_frac{f} = cell( Fracture(g).N_Total , 1 );
                        
                        Fracture(f).     Af_frac_EDFM{g} = cell( Fracture(f).N_Total , 1 );
                        Fracture(f).avgDist_frac_EDFM{g} = cell( Fracture(f).N_Total , 1 );
                        Fracture(f).     CI_frac_EDFM{g} = cell( Fracture(f).N_Total , 1 );
                        Fracture(f). Af_frac_EDFM_sum{g} = 0;
                        
                        Fracture(g).     Af_frac_EDFM{f} = cell( Fracture(g).N_Total , 1 );
                        Fracture(g).avgDist_frac_EDFM{f} = cell( Fracture(g).N_Total , 1 );
                        Fracture(g).     CI_frac_EDFM{f} = cell( Fracture(g).N_Total , 1 );
                        Fracture(g). Af_frac_EDFM_sum{f} = 0;
                        
                        % Loop for fracture f
                        If_List = find(Fracture(f).Overlap_frac{g});
                        for nf = 1:length(If_List)
                            If = If_List(nf);
                            [i_f,j_f,~] = Index_Global_To_Local_3D(Fracture(f).N_Length_AB, Fracture(f).N_Width_AD, 1, If);
                            
                            % Corner Points of the 1st fracture cell
                            PointA = [ Fracture(f).GridCoords(i_f  ,j_f  ,1) ; Fracture(f).GridCoords(i_f  ,j_f  ,2) ; Fracture(f).GridCoords(i_f  ,j_f  ,3) ];
                            PointB = [ Fracture(f).GridCoords(i_f  ,j_f+1,1) ; Fracture(f).GridCoords(i_f  ,j_f+1,2) ; Fracture(f).GridCoords(i_f  ,j_f+1,3) ];
                            PointC = [ Fracture(f).GridCoords(i_f+1,j_f+1,1) ; Fracture(f).GridCoords(i_f+1,j_f+1,2) ; Fracture(f).GridCoords(i_f+1,j_f+1,3) ];
                            PointD = [ Fracture(f).GridCoords(i_f+1,j_f  ,1) ; Fracture(f).GridCoords(i_f+1,j_f  ,2) ; Fracture(f).GridCoords(i_f+1,j_f  ,3) ];
                            Plane_fracCell_1 = planeSegment_FracGen(PointA,PointB,PointC,PointD);
                            
                            % Loop for fracture g
                            Ind_f_Intersect  = 0;
                            Ig_List = find(Fracture(g).Overlap_frac{f});
                            for ng = 1:length(Ig_List)
                                Ig = Ig_List(ng);
                                [i_g,j_g,~] = Index_Global_To_Local_3D(Fracture(g).N_Length_AB, Fracture(g).N_Width_AD, 1, Ig);
                                
                                % Corner Points of the 2nd fracture cell
                                PointA = [ Fracture(g).GridCoords(i_g  ,j_g  ,1) ; Fracture(g).GridCoords(i_g  ,j_g  ,2) ; Fracture(g).GridCoords(i_g  ,j_g  ,3) ];
                                PointB = [ Fracture(g).GridCoords(i_g  ,j_g+1,1) ; Fracture(g).GridCoords(i_g  ,j_g+1,2) ; Fracture(g).GridCoords(i_g  ,j_g+1,3) ];
                                PointC = [ Fracture(g).GridCoords(i_g+1,j_g+1,1) ; Fracture(g).GridCoords(i_g+1,j_g+1,2) ; Fracture(g).GridCoords(i_g+1,j_g+1,3) ];
                                PointD = [ Fracture(g).GridCoords(i_g+1,j_g  ,1) ; Fracture(g).GridCoords(i_g+1,j_g  ,2) ; Fracture(g).GridCoords(i_g+1,j_g  ,3) ];
                                Plane_fracCell_2 = planeSegment_FracGen(PointA,PointB,PointC,PointD);
                                
                                [~, intersectPoints] = Plane_fracCell_1.Obtain_PlaneSegment_PlaneSegment_Intersection( Plane_fracCell_2, Epsilon );
                                
                                intersectCoord_frac_final = [];
                                
                                if ~isempty(intersectPoints)
                                    for nr = 1 : size(intersectPoints,2)
                                        
                                        % Removing the Points that are not inside the fracture rectangle
                                        isInside1 = Plane_fracCell_1.Is_Point_Inside_PlaneSegment(intersectPoints(:,nr), Epsilon);
                                        isInside2 = Plane_fracCell_2.Is_Point_Inside_PlaneSegment(intersectPoints(:,nr), Epsilon);
                                        if ( isInside1 == 0 ) || ( isInside2 == 0 )
                                            intersectPoints(:,nr) = [ NaN ; NaN ; NaN ];
                                            continue;
                                        end
                                        
                                        % Removing the Points that are repeated
                                        subtracted = ( intersectPoints - intersectPoints(:,nr) );
                                        subtracted = sqrt( subtracted(1,:).^2 + subtracted(2,:).^2 + subtracted(3,:).^2 );
                                        if ~isempty( find( subtracted([1:nr-1 , nr+1:end]) < Epsilon , 1 ) )
                                            intersectPoints(:,nr) = [ NaN ; NaN ; NaN ];
                                            continue;
                                        end
                                        
                                        if ~isnan(intersectPoints(1,nr))
                                            intersectCoord_frac_final = [ intersectCoord_frac_final , intersectPoints(:,nr) ];
                                        end
                                    end
                                    
                                    if ( size( intersectCoord_frac_final , 2 ) > 2 )
                                        error('More than two intersection Points exist between two intersecting fracture cells! The implementation is wrong!');
                                    end
                                    
                                    % Obtaining the area fractions and average distances and writing them into structure
                                    if ( size( intersectCoord_frac_final , 2 ) == 2 )
                                        Ind_f_Intersect = Ind_f_Intersect + 1;
                                        Ind_g_Intersect = size(Fracture(g).intersectCoord_frac{f}{Ig},1) + 1;
                                        Fracture(f).intersectCoord_frac{g}{If}{Ind_f_Intersect,1} = Ig;
                                        Fracture(f).intersectCoord_frac{g}{If}{Ind_f_Intersect,2} = intersectCoord_frac_final;
                                        Fracture(g).intersectCoord_frac{f}{Ig}{Ind_g_Intersect,1} = If;
                                        Fracture(g).intersectCoord_frac{f}{Ig}{Ind_g_Intersect,2} = intersectCoord_frac_final;
                                        
                                        Fracture(f).NumOf_fracConn_EDFM(If) = Fracture(f).NumOf_fracConn_EDFM(If) + 1;
                                        Fracture(g).NumOf_fracConn_EDFM(Ig) = Fracture(g).NumOf_fracConn_EDFM(Ig) + 1;
                                        
                                        Fracture(f).     Af_frac_EDFM{g}{If}(Ind_f_Intersect,1) = Ig;
                                        Fracture(f).avgDist_frac_EDFM{g}{If}(Ind_f_Intersect,1) = Ig;
                                        Fracture(f).     CI_frac_EDFM{g}{If}(Ind_f_Intersect,1) = Ig;
                                        Fracture(g).     Af_frac_EDFM{f}{Ig}(Ind_g_Intersect,1) = If;
                                        Fracture(g).avgDist_frac_EDFM{f}{Ig}(Ind_g_Intersect,1) = If;
                                        Fracture(g).     CI_frac_EDFM{f}{Ig}(Ind_g_Intersect,1) = If;
                                        
                                        intersectLine = lineSegment_FracGen( intersectCoord_frac_final(:,1) , intersectCoord_frac_final(:,2) );
                                        intersectLine.AddPointM;
                                        
                                        [Af_frac1 , avgDist_frac1, Collinearity1 ] = obj.Line_Plane_Connectivity_3D( Plane_fracCell_1 , intersectLine );
                                        [Af_frac2 , avgDist_frac2, Collinearity2 ] = obj.Line_Plane_Connectivity_3D( Plane_fracCell_2 , intersectLine );

                                        if ( Collinearity1 == 1 ),  Af_frac2 = Af_frac2 /2;  end
                                        if ( Collinearity2 == 1 ),  Af_frac1 = Af_frac1 /2;  end
                                        
                                        Fracture(f).     Af_frac_EDFM{g}{If}(Ind_f_Intersect,2) = Af_frac1;
                                        Fracture(f).avgDist_frac_EDFM{g}{If}(Ind_f_Intersect,2) = avgDist_frac1;
                                        Fracture(g).     Af_frac_EDFM{f}{Ig}(Ind_g_Intersect,2) = Af_frac2;
                                        Fracture(g).avgDist_frac_EDFM{f}{Ig}(Ind_g_Intersect,2) = avgDist_frac2;
                                        
                                        % Each two none-neighboring cells
                                        % have separate calculation for
                                        % connectivity index (CI). However,
                                        % they need to share one unique CI.
                                        % Therefore both CIs will be first
                                        % harmonically averaged and then
                                        % will be added to the dataset.
                                        CI_frac1 = Af_frac1 / avgDist_frac1;
                                        CI_frac2 = Af_frac2 / avgDist_frac2;
                                        CI_Harmonic = 2 * CI_frac1 * CI_frac2 / (CI_frac1 + CI_frac2);
                                        
                                        Fracture(f).CI_frac_EDFM{g}{If}(Ind_f_Intersect,2) = CI_Harmonic;
                                        Fracture(g).CI_frac_EDFM{f}{Ig}(Ind_g_Intersect,2) = CI_Harmonic;
                                        
                                        Fracture(f).Af_frac_EDFM_sum{g} = Fracture(f).Af_frac_EDFM_sum{g} + Af_frac1;
                                        Fracture(g).Af_frac_EDFM_sum{f} = Fracture(g).Af_frac_EDFM_sum{f} + Af_frac2;
                                    end
                                end
                            end
                        end
                    end
                end
            end
            fprintf('\n---------------------------------------------------------\n');
        end
        function Add_pEDFM_Connectivities(obj)
            Fracture = obj.FracturesGrid.Fracture;
            Epsilon = obj.ReservoirGrid.Epsilon;
            %% Assigning reservoir properties
            % Number of reservoir grid cells in each direction
            NX = obj.ReservoirGrid.NX;
            NY = obj.ReservoirGrid.NY;
            NZ = obj.ReservoirGrid.NZ;
            % Coordinates of reservoir cell centers in each direction
            Xcm = obj.ReservoirGrid.Xcm;
            Ycm = obj.ReservoirGrid.Ycm;
            Zcm = obj.ReservoirGrid.Zcm;
            % Coordinates of reservoir grid nodes in each direction
            Xim = obj.ReservoirGrid.Xim;
            Yim = obj.ReservoirGrid.Yim;
            Zim = obj.ReservoirGrid.Zim;
            %% Adding pEDFM Connectivities
            fprintf('Adding pEDFM Connectivities:\n');
            fprintf('---> Fracture ');
            
            % Looping through fractures
            for f = 1 : length(Fracture)
                if (f>1),  fprintf(repmat('\b', 1, 5+27));  end
                fprintf('%02d/%02d',f,length(Fracture));
                fprintf(' ---> Grid cell ');
                for If = 1 : Fracture(f).N_Total
                    if (If>1),  fprintf(repmat('\b', 1, 11));  end
                    fprintf('%05d/%05d',If,Fracture(f).N_Total);
                    
                    %% 1- Obtaining pEDFM connectivities of this frature with matrix
                    Af = sortrows(Fracture(f).Af_rock_EDFM{If});
                    Im_List = Af(:,1);
                    % Initializing the necessary variable to collect the conectivities data
                    counter = 0;
                    pEDFM_List = [];
                    pEDFM_Af = [];
                    pEDFM_AvgDist = [];
                    
                    % Looping over the overlapped matrix cells
                    for m = 1 : length(Im_List)
                        Af_rock_EDFM_New = 0;
                        Im = Im_List(m);
                        [i,j,k] = Index_Global_To_Local_3D(NX, NY, NZ, Im);
                        
                        % Finding the possible intersection point between
                        % the fracture plane and the 3 axis lines the pass
                        % throught the cell center of this matrix cell.
                        CellCenter_This = [Xcm(i);Ycm(j);Zcm(k)];
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % Checking left and right
                        LineAlongX = lineInf_FracGen( [1;0;0] );
                        LineAlongX.AddEquation(CellCenter_This);
                        [Geostatus, intersectPoint] = Fracture(f).Obtain_LineInf_PlaneInf_Intersection(LineAlongX, Epsilon);
                        
                        if Geostatus.haveIntersect == 1 && ~isempty(intersectPoint) && NX>1
                            % Finding the cell(s) involved in this intersection
                            [minVal,Ind_Center] = min(abs(Xcm - intersectPoint(1)));
                            if strcmp(Fracture(f).isParallelToRefSurface,'AlongYZ') && (minVal < Epsilon) && (Ind_Center == i)
                                i_center = [i-1 ; i+1];
                                i_interface = [i ; i+1];
                                % If the fracture is parallel AlongYZ and cuts throught the cell center,
                                % the connectivities for both left and right neighbors should be changed.
                            else
                                [~,i_interface] = min(abs(Xim - intersectPoint(1)));
                                if i_interface > 1 && i_interface < NX+1
                                    if i_interface > i
                                        i_center = i_interface;
                                    else
                                        i_center = i_interface - 1;
                                    end
                                end
                            end
                            for ii = 1 : 1%length(i_center)
                                if i_interface(ii) > 1 && i_interface(ii) < NX+1
                                    Im_intersected = Index_Local_To_Global_3D(NX, NY, NZ, i_center(ii), j, k);
                                    
                                    % Check if the intersected cell is already on the list or not
                                    if ~ismember(Im_intersected,pEDFM_List)
                                        pEDFM_List = [pEDFM_List;Im_intersected];
                                        pEDFM_Af = [pEDFM_Af;0];
                                        pEDFM_AvgDist = [pEDFM_AvgDist;0];
                                        counter = counter+1;
                                        Index = counter;
                                    else
                                        Index = find(pEDFM_List==Im_intersected);
                                    end
                                    
                                    % Calculating the projected area fraction Af_proj:
                                    % the Cosine of the projection angle can be
                                    % calculated as dot product between the
                                    % normal vectors of fracture plate and
                                    % projection face.
                                    Af_proj = (Af(m,2) / 2) * abs(dot(Fracture(f).n_vec,[1;0;0]));
                                    pEDFM_Af(Index) = pEDFM_Af(Index) + Af_proj;
                                    
                                    % The Af_rock_EDFM also needs to change
                                    Af_rock_EDFM_New = Af_rock_EDFM_New + Af_proj;
                                    
                                    % Calculating average distance AvgDist_proj:
                                    CubeCornerStart = [ Xim(i_center(ii))   ; Yim(j)   ; Zim(k)   ];
                                    CubeCornerEnd   = [ Xim(i_center(ii)+1) ; Yim(j+1) ; Zim(k+1) ];
                                    AvgDist_proj = Average_Distance_Cube_From_Plane( Fracture(f).Equation, CubeCornerStart , CubeCornerEnd , 5 );
                                    pEDFM_AvgDist(Index) = AvgDist_proj;
                                    
                                    % Measuring the alpha:
                                    % alpha is a the fraction of blockage (from 0 to 1) between
                                    % neighboring matrix cells due to existence of fracture plates when
                                    % considering pEDFM method. the alpha=0 means no blockage and
                                    % alpha=1 means total blockage.
                                    alpha = (Af_proj) / (obj.ReservoirGrid.DY * obj.ReservoirGrid.DZ);
                                    %if (alpha>1) || (alpha<0), error('alpha cannot be bigger than 1'); end
                                    obj.ReservoirGrid.alpha_Tx(i_interface(ii),j,k) = obj.ReservoirGrid.alpha_Tx(i_interface(ii),j,k) + alpha;
                                    obj.ReservoirGrid.alpha_Tx(i_interface(ii),j,k) = min( obj.ReservoirGrid.alpha_Tx(i_interface(ii),j,k) , 1 );
                                end
                            end
                        end
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % Checking back and front
                        LineAlongY = lineInf_FracGen( [0;1;0] );
                        LineAlongY.AddEquation(CellCenter_This);
                        [Geostatus, intersectPoint] = Fracture(f).Obtain_LineInf_PlaneInf_Intersection(LineAlongY, Epsilon);
                        
                        if Geostatus.haveIntersect == 1 && ~isempty(intersectPoint) && NY>1
                            % Finding the cell(s) involved in this intersection
                            [minVal,Ind_Center] = min(abs(Ycm - intersectPoint(2)));
                            if strcmp(Fracture(f).isParallelToRefSurface,'AlongXZ') && (minVal < Epsilon) && (Ind_Center == j)
                                j_center = [j-1 ; j+1];
                                j_interface = [j ; j+1];
                                % If the fracture is parallel AlongXZ and cuts throught the cell center,
                                % the connectivities for both back and front neighbors should be changed.
                            else
                                [~,j_interface] = min(abs(Yim - intersectPoint(2)));
                                if j_interface > 1 && j_interface < NY+1
                                    if j_interface > j
                                        j_center = j_interface;
                                    else
                                        j_center = j_interface - 1;
                                    end
                                end
                            end
                            for jj = 1 : 1%length(j_center)
                                if j_interface(jj) > 1 && j_interface(jj) < NY+1
                                    Im_intersected = Index_Local_To_Global_3D(NX, NY, NZ, i, j_center(jj), k);
                                    
                                    % Check if the intersected cell is already on the list or not
                                    if ~ismember(Im_intersected,pEDFM_List)
                                        pEDFM_List = [pEDFM_List;Im_intersected];
                                        pEDFM_Af = [pEDFM_Af;0];
                                        pEDFM_AvgDist = [pEDFM_AvgDist;0];
                                        counter = counter+1;
                                        Index = counter;
                                    else
                                        Index = find(pEDFM_List==Im_intersected);
                                    end
                                    
                                    % Calculating the projected area fraction Af_proj:
                                    % the Cosine of the projection angle can be
                                    % calculated as dot product between the
                                    % normal vectors of fracture plate and
                                    % projection face.
                                    Af_proj = (Af(m,2) / 2) * abs(dot(Fracture(f).n_vec,[0;1;0]));
                                    pEDFM_Af(Index) = pEDFM_Af(Index) + Af_proj;
                                    
                                    % The Af_rock_EDFM also needs to change
                                    Af_rock_EDFM_New = Af_rock_EDFM_New + Af_proj;
                                    
                                    % Calculating average distance AvgDist_proj:
                                    CubeCornerStart = [ Xim(i)   ; Yim(j_center(jj))   ; Zim(k)   ];
                                    CubeCornerEnd   = [ Xim(i+1) ; Yim(j_center(jj)+1) ; Zim(k+1) ];
                                    AvgDist_proj = Average_Distance_Cube_From_Plane( Fracture(f).Equation, CubeCornerStart , CubeCornerEnd , 5 );
                                    pEDFM_AvgDist(Index) = AvgDist_proj;
                                    
                                    % Measuring the alpha:
                                    % alpha is a the fraction of blockage (from 0 to 1) between
                                    % neighboring matrix cells due to existence of fracture plates when
                                    % considering pEDFM method. the alpha=0 means no blockage and
                                    % alpha=1 means total blockage.
                                    alpha = (Af_proj) / (obj.ReservoirGrid.DX * obj.ReservoirGrid.DZ);
                                    %if (alpha>1) || (alpha<0), error('alpha cannot be bigger than 1'); end
                                    obj.ReservoirGrid.alpha_Ty(i,j_interface(jj),k) = obj.ReservoirGrid.alpha_Ty(i,j_interface(jj),k) + alpha;
                                    obj.ReservoirGrid.alpha_Ty(i,j_interface(jj),k) = min( obj.ReservoirGrid.alpha_Ty(i,j_interface(jj),k) , 1 );
                                end
                            end
                        end
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % Checking bottom and top
                        LineAlongZ = lineInf_FracGen( [0;0;1] );
                        LineAlongZ.AddEquation(CellCenter_This);
                        [Geostatus, intersectPoint] = Fracture(f).Obtain_LineInf_PlaneInf_Intersection(LineAlongZ, Epsilon);
                        
                        if Geostatus.haveIntersect == 1 && ~isempty(intersectPoint) && NZ>1
                            % Finding the cell(s) involved in this intersection
                            [minVal,Ind_Center] = min(abs(Zcm - intersectPoint(3)));
                            if strcmp(Fracture(f).isParallelToRefSurface,'AlongXY') && (minVal < Epsilon) && (Ind_Center == k)
                                k_center = [k-1 ; k+1];
                                k_interface = [k ; k+1];
                                % If the fracture is parallel AlongXY and cuts throught the cell center,
                                % the connectivities for both bottom and top neighbors should be changed.
                            else
                                [~,k_interface] = min(abs(Zim - intersectPoint(3)));
                                if k_interface > 1 && k_interface < NZ+1
                                    if k_interface > k
                                        k_center = k_interface;
                                    else
                                        k_center = k_interface - 1;
                                    end
                                end
                            end
                            for kk = 1 : 1%length(k_center)
                                if k_interface(kk) > 1 && k_interface(kk) < NZ+1
                                    Im_intersected = Index_Local_To_Global_3D(NX, NY, NZ, i, j, k_center(kk));
                                    
                                    % Check if the intersected cell is already on the list or not
                                    if ~ismember(Im_intersected,pEDFM_List)
                                        pEDFM_List = [pEDFM_List;Im_intersected];
                                        pEDFM_Af = [pEDFM_Af;0];
                                        pEDFM_AvgDist = [pEDFM_AvgDist;0];
                                        counter = counter+1;
                                        Index = counter;
                                    else
                                        Index = find(pEDFM_List==Im_intersected);
                                    end
                                    
                                    % Calculating the projected area fraction Af_proj:
                                    % the Cosine of the projection angle can be
                                    % calculated as dot product between the
                                    % normal vectors of fracture plate and
                                    % projection face.
                                    Af_proj = (Af(m,2) / 2) * abs(dot(Fracture(f).n_vec,[0;0;1]));
                                    pEDFM_Af(Index) = pEDFM_Af(Index) + Af_proj;
                                    
                                    % The Af_rock_EDFM also needs to change
                                    Af_rock_EDFM_New = Af_rock_EDFM_New + Af_proj;
                                    
                                    % Calculating average distance AvgDist_proj:
                                    CubeCornerStart = [ Xim(i)   ; Yim(j)   ; Zim(k_center(kk))   ];
                                    CubeCornerEnd   = [ Xim(i+1) ; Yim(j+1) ; Zim(k_center(kk)+1) ];
                                    AvgDist_proj = Average_Distance_Cube_From_Plane( Fracture(f).Equation, CubeCornerStart , CubeCornerEnd , 5 );
                                    pEDFM_AvgDist(Index) = AvgDist_proj;
                                    
                                    % Measuring the alpha:
                                    % alpha is a the fraction of blockage (from 0 to 1) between
                                    % neighboring matrix cells due to existence of fracture plates when
                                    % considering pEDFM method. the alpha=0 means no blockage and
                                    % alpha=1 means total blockage.
                                    alpha = (Af_proj) / (obj.ReservoirGrid.DX * obj.ReservoirGrid.DY);
                                    %if (alpha>1) || (alpha<0), error('alpha cannot be bigger than 1'); end
                                    obj.ReservoirGrid.alpha_Tz(i,j,k_interface(kk)) = obj.ReservoirGrid.alpha_Tz(i,j,k_interface(kk)) + alpha;
                                    obj.ReservoirGrid.alpha_Tz(i,j,k_interface(kk)) = min( obj.ReservoirGrid.alpha_Tz(i,j,k_interface(kk)) , 1 );
                                end
                            end
                        end
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % Modifying the Af and CI of EDFM connectivities of
                        % cell "Im" due to pEDFM connectivities of neighboring cells
                        if Af_rock_EDFM_New ~= 0
                            ind = find( Fracture(f).Af_rock_EDFM{If}(:,1) == Im );
                            if Fracture(f).Af_rock_EDFM{If}(ind,2) < Af_rock_EDFM_New
                                warning('Fracture(%2d).Af_rock_EDFM{%2d}(%2d,2) is smaller than Af_rock_EDFM_New!',f,If,ind);
                            end
                            Fracture(f).Af_rock_EDFM{If}(ind,2) = Af_rock_EDFM_New;
                            Fracture(f).CI_rock_EDFM{If}(ind,2) = Fracture(f).Af_rock_EDFM{If}(ind,2) / Fracture(f).avgDist_rock_EDFM{If}(ind,2);
                        end
                    end
                    % Adding the pEDFM information to the fracture
                    pEDFM_Temp = [pEDFM_List,pEDFM_Af,pEDFM_AvgDist,pEDFM_Af./pEDFM_AvgDist];
                    pEDFM_Temp = sortrows(pEDFM_Temp,1);
                    Fracture(f).     Af_rock_pEDFM{If} = pEDFM_Temp(:,[1,2]);
                    Fracture(f).avgDist_rock_pEDFM{If} = pEDFM_Temp(:,[1,3]);
                    Fracture(f).     CI_rock_pEDFM{If} = pEDFM_Temp(:,[1,4]);
                    
                    %% 2- Obtaining pEDFM connectivities of this frature with other fractures
                    if Fracture(f).NumOf_fracConn_EDFM(If) > 0
                        for g = 1 : length(Fracture)
                            if g==f,  continue;  end
                            if isempty(Fracture(f).Af_frac_EDFM{g}),  continue;  end
                            if ~isempty(Fracture(f).Af_frac_EDFM{g}{If})
                                if isempty(Fracture(f).Af_frac_pEDFM{g})
                                    Fracture(f).     Af_frac_pEDFM{g} = cell(Fracture(f).N_Total,1);
                                    Fracture(f).avgDist_frac_pEDFM{g} = cell(Fracture(f).N_Total,1);
                                    Fracture(f).     CI_frac_pEDFM{g} = cell(Fracture(f).N_Total,1);
                                end
                                Aff = sortrows(Fracture(f).Af_frac_EDFM{g}{If});
                                Ig_List = Aff(:,1);
                                counter = 0;
                                pEDFM_List = [];
                                pEDFM_Af = [];
                                pEDFM_AvgDist = [];
                                IntersectLine = lineSegment_FracGen( Fracture(f).intersectCoord_fracturePlate{g}(:,2),Fracture(f).intersectCoord_fracturePlate{g}(:,1) );
                                IntersectLine.AddPointM;
                                for n = 1 : length(Ig_List)
                                    Ig = Ig_List(n);
                                    [i_g,j_g,~] = Index_Global_To_Local_3D(Fracture(g).N_Length_AB, Fracture(g).N_Width_AD  , 1, Ig);
                                    CellCenter_This = Fracture(g).CellCenterCoordsV2(Ig,:)';
                                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                    % Checking along length
                                    LineAlongLength = lineInf_FracGen( Fracture(g).AB_vec );
                                    LineAlongLength.AddPoint0( CellCenter_This );
                                    [Geostatus, intersectPoint] = IntersectLine.Obtain_LineInf_LineInf_Intersection( LineAlongLength, Epsilon);
                                    
                                    if Geostatus.haveIntersect == 1 && ~isempty(intersectPoint) && Fracture(g).N_Length_AB>1
                                        % Finding the cell(s) involved in this intersection
                                        CellCenterCoords = zeros(Fracture(g).N_Length_AB,3);
                                        CellCenterCoords(:,1) = Fracture(g).CellCenterCoordsV1(:,j_g,1);
                                        CellCenterCoords(:,2) = Fracture(g).CellCenterCoordsV1(:,j_g,2);
                                        CellCenterCoords(:,3) = Fracture(g).CellCenterCoordsV1(:,j_g,3);
                                        CellCenterDistance = CellCenterCoords - repmat(intersectPoint',Fracture(g).N_Length_AB,1);
                                        [minVal,Ind_Center] = min(vecnorm(CellCenterDistance'));
                                        if ( norm( dot(IntersectLine.unit_vec, LineAlongLength.unit_vec) ) < Epsilon )&& (minVal < Epsilon) && (Ind_Center == i_g)
                                            i_g_center = [i_g-1 ; i_g+1];
                                            i_g_interface = [i_g ; i_g+1];
                                            % If the intersection line is parallel to length of fracture(g) and cuts throught the cell center,
                                            % the connectivities for both neighbors should be changed.
                                        else
                                            InterfaceCoords = CellCenterCoords - repmat([Fracture(g).D_Length_AB_vec/2]',Fracture(g).N_Length_AB,1);
                                            InterfaceCoords(end+1 , :) = InterfaceCoords(end,:)+Fracture(g).D_Length_AB_vec';
                                            InterfaceDistance = InterfaceCoords - repmat(intersectPoint',Fracture(g).N_Length_AB+1,1);
                                            [~,i_g_interface] = min(vecnorm(InterfaceDistance'));
                                            if i_g_interface > 1 && i_g_interface < Fracture(g).N_Length_AB+1
                                                if i_g_interface > i_g
                                                    i_g_center = i_g_interface;
                                                else
                                                    i_g_center = i_g_interface - 1;
                                                end
                                            end
                                        end
                                        for ii_g = 1 : 1 %length(i_g_center)
                                            if i_g_interface(ii_g) > 1 && i_g_interface(ii_g) < Fracture(g).N_Length_AB+1
                                                Ig_intersected = Index_Local_To_Global_3D(Fracture(g).N_Length_AB, Fracture(g).N_Width_AD, 1, i_g_center(ii_g), j_g, 1);
                                                
                                                % Check if the intersected cell is already on the list or not
                                                if ~ismember(Ig_intersected,pEDFM_List)
                                                    pEDFM_List = [pEDFM_List;Ig_intersected];
                                                    pEDFM_Af = [pEDFM_Af;0];
                                                    pEDFM_AvgDist = [pEDFM_AvgDist;0];
                                                    counter = counter+1;
                                                    Index = counter;
                                                else
                                                    Index = find(pEDFM_List==Ig_intersected);
                                                end
                                                % Calculating the projected area fraction Aff_proj:
                                                % the Cosine of the projection angle can be
                                                % calculated as dot product between the
                                                % unit vectors of fracture length and
                                                % projection line.
                                                Aff_proj = (Aff(n,2) / 2) * abs(dot( Fracture(g).AD_vec/norm(Fracture(g).AD_vec) , IntersectLine.unit_vec) );
                                                pEDFM_Af(Index) = pEDFM_Af(Index) + Aff_proj;
                                                
                                                % Calculating average distance AvgDist_proj:
                                                PointA = [ Fracture(g).GridCoords(i_g_center(ii_g)  ,j_g+1,1)
                                                           Fracture(g).GridCoords(i_g_center(ii_g)  ,j_g+1,2)
                                                           Fracture(g).GridCoords(i_g_center(ii_g)  ,j_g+1,3) ];
                                                PointB = [ Fracture(g).GridCoords(i_g_center(ii_g)  ,j_g  ,1)
                                                           Fracture(g).GridCoords(i_g_center(ii_g)  ,j_g  ,2)
                                                           Fracture(g).GridCoords(i_g_center(ii_g)  ,j_g  ,3) ];
                                                PointC = [ Fracture(g).GridCoords(i_g_center(ii_g)+1,j_g  ,1)
                                                           Fracture(g).GridCoords(i_g_center(ii_g)+1,j_g  ,2)
                                                           Fracture(g).GridCoords(i_g_center(ii_g)+1,j_g  ,3) ];
                                                PointD = [ Fracture(g).GridCoords(i_g_center(ii_g)+1,j_g+1,1)
                                                           Fracture(g).GridCoords(i_g_center(ii_g)+1,j_g+1,2)
                                                           Fracture(g).GridCoords(i_g_center(ii_g)+1,j_g+1,3) ];
                                                NeighborCell = planeSegment_FracGen(PointA,PointB,PointC,PointD);
                                                AvgDist_proj = NeighborCell.Obtain_Average_Dsitance_LineInf_From_PlaneSegment(IntersectLine, 6);
                                                pEDFM_AvgDist(Index) = AvgDist_proj;
                                                
                                                % Measuring the alpha:
                                                % alpha is a the fraction of blockage (from 0 to 1) between
                                                % neighboring matrix cells due to existence of fracture plates when
                                                % considering pEDFM method. the alpha=0 means no blockage and
                                                % alpha=1 means total blockage.
                                                alpha = (Aff_proj) / Fracture(g).D_Width_AD;
                                                Fracture(g).alpha_Tx(i_g_interface(ii_g),j_g) = Fracture(g).alpha_Tx(i_g_interface(ii_g),j_g) + alpha;
                                                Fracture(g).alpha_Tx(i_g_interface(ii_g),j_g) = min( Fracture(g).alpha_Tx(i_g_interface(ii_g),j_g) , 1 );
                                            end
                                        end
                                    end
                                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                    % Checking along width
                                    LineAlongWidth = lineInf_FracGen( Fracture(g).AD_vec );
                                    LineAlongWidth.AddPoint0( CellCenter_This );
                                    [Geostatus, intersectPoint] = IntersectLine.Obtain_LineInf_LineInf_Intersection( LineAlongWidth, Epsilon);
                                    
                                    if Geostatus.haveIntersect == 1 && ~isempty(intersectPoint) && Fracture(g).N_Width_AD>1
                                        % Finding the cell(s) involved in this intersection
                                        CellCenterCoords = zeros(Fracture(g).N_Width_AD,3);
                                        CellCenterCoords(:,1) = Fracture(g).CellCenterCoordsV1(i_g,:,1)';
                                        CellCenterCoords(:,2) = Fracture(g).CellCenterCoordsV1(i_g,:,2)';
                                        CellCenterCoords(:,3) = Fracture(g).CellCenterCoordsV1(i_g,:,3)';
                                        CellCenterDistance = CellCenterCoords - repmat(intersectPoint',Fracture(g).N_Width_AD,1);
                                        [minVal,Ind_Center] = min(vecnorm(CellCenterDistance'));
                                        if ( norm( dot(IntersectLine.unit_vec, LineAlongWidth.unit_vec) ) < Epsilon )&& (minVal < Epsilon) && (Ind_Center == j_g)
                                            j_g_center = [j_g-1 ; j_g+1];
                                            j_g_interface = [j_g ; j_g+1];
                                            % If the intersection line is parallel to width of fracture(g) and cuts throught the cell center,
                                            % the connectivities for both neighbors should be changed.
                                        else
                                            InterfaceCoords = CellCenterCoords - repmat([Fracture(g).D_Width_AD_vec/2]',Fracture(g).N_Width_AD,1);
                                            InterfaceCoords(end+1 , :) = InterfaceCoords(end,:)+Fracture(g).D_Width_AD_vec';
                                            InterfaceDistance = InterfaceCoords - repmat(intersectPoint',Fracture(g).N_Width_AD+1,1);
                                            [~,j_g_interface] = min(vecnorm(InterfaceDistance'));
                                            if j_g_interface > 1 && j_g_interface < Fracture(g).N_Width_AD+1
                                                if j_g_interface > j_g
                                                    j_g_center = j_g_interface;
                                                else
                                                    j_g_center = j_g_interface - 1;
                                                end
                                            end
                                        end
                                        for jj_g = 1 : 1 %length(j_g_center)
                                            if j_g_interface(jj_g) > 1 && j_g_interface(jj_g) < Fracture(g).N_Width_AD+1
                                                Ig_intersected = Index_Local_To_Global_3D(Fracture(g).N_Length_AB, Fracture(g).N_Width_AD, 1, i_g, j_g_center(jj_g), 1);
                                                
                                                % Check if the intersected cell is already on the list or not
                                                if ~ismember(Ig_intersected,pEDFM_List)
                                                    pEDFM_List = [pEDFM_List;Ig_intersected];
                                                    pEDFM_Af = [pEDFM_Af;0];
                                                    pEDFM_AvgDist = [pEDFM_AvgDist;0];
                                                    counter = counter+1;
                                                    Index = counter;
                                                else
                                                    Index = find(pEDFM_List==Ig_intersected);
                                                end
                                                % Calculating the projected area fraction Aff_proj:
                                                % the Cosine of the projection angle can be
                                                % calculated as dot product between the
                                                % unit vectors of fracture width and
                                                % projection line.
                                                Aff_proj = (Aff(n,2) / 2) * abs(dot( Fracture(g).AB_vec/norm(Fracture(g).AB_vec) , IntersectLine.unit_vec) );
                                                pEDFM_Af(Index) = pEDFM_Af(Index) + Aff_proj;
                                                
                                                % Calculating average distance AvgDist_proj:
                                                PointA = [ Fracture(g).GridCoords(i_g  ,j_g_center(jj_g)+1,1)
                                                           Fracture(g).GridCoords(i_g  ,j_g_center(jj_g)+1,2)
                                                           Fracture(g).GridCoords(i_g  ,j_g_center(jj_g)+1,3) ];
                                                PointB = [ Fracture(g).GridCoords(i_g  ,j_g_center(jj_g)  ,1)
                                                           Fracture(g).GridCoords(i_g  ,j_g_center(jj_g)  ,2)
                                                           Fracture(g).GridCoords(i_g  ,j_g_center(jj_g)  ,3) ];
                                                PointC = [ Fracture(g).GridCoords(i_g+1,j_g_center(jj_g)  ,1)
                                                           Fracture(g).GridCoords(i_g+1,j_g_center(jj_g)  ,2)
                                                           Fracture(g).GridCoords(i_g+1,j_g_center(jj_g)  ,3) ];
                                                PointD = [ Fracture(g).GridCoords(i_g+1,j_g_center(jj_g)+1,1)
                                                           Fracture(g).GridCoords(i_g+1,j_g_center(jj_g)+1,2)
                                                           Fracture(g).GridCoords(i_g+1,j_g_center(jj_g)+1,3) ];
                                                NeighborCell = planeSegment_FracGen(PointA,PointB,PointC,PointD);
                                                AvgDist_proj = NeighborCell.Obtain_Average_Dsitance_LineInf_From_PlaneSegment(IntersectLine, 6);
                                                pEDFM_AvgDist(Index) = AvgDist_proj;
                                                
                                                % Measuring the alpha:
                                                % alpha is a the fraction of blockage (from 0 to 1) between
                                                % neighboring matrix cells due to existence of fracture plates when
                                                % considering pEDFM method. the alpha=0 means no blockage and
                                                % alpha=1 means total blockage.
                                                alpha = (Aff_proj) / Fracture(g).D_Length_AB;
                                                Fracture(g).alpha_Ty(i_g,j_g_interface(jj_g)) = Fracture(g).alpha_Ty(i_g,j_g_interface(jj_g)) + alpha;
                                                Fracture(g).alpha_Ty(i_g,j_g_interface(jj_g)) = min( Fracture(g).alpha_Ty(i_g,j_g_interface(jj_g)) , 1 );
                                            end
                                        end
                                    end
                                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                end
                                % Adding the pEDFM information to the fracture
                                pEDFM_Temp = [pEDFM_List,pEDFM_Af,pEDFM_AvgDist,pEDFM_Af./pEDFM_AvgDist];
                                if ~isempty(pEDFM_Temp)
                                    pEDFM_Temp = sortrows(pEDFM_Temp,1);
                                    Fracture(f).Af_frac_pEDFM{g}{If}      = pEDFM_Temp(:,[1,2]);
                                    Fracture(f).avgDist_frac_pEDFM{g}{If} = pEDFM_Temp(:,[1,3]);
                                    Fracture(f).CI_frac_pEDFM{g}{If}      = pEDFM_Temp(:,[1,4]);
                                    Fracture(f).NumOf_fracConn_pEDFM(If)  = Fracture(f).NumOf_fracConn_pEDFM(If) + size(pEDFM_Temp,1);
                                end
                            end
                        end
                    end
                end
            end
            
            % Each two none-neighboring cells have separate calculation for
            % connectivity index (CI). However, they need to share one 
            % unique CI. Therefore both CIs (if they co-exist) will be
            % first harmonically averaged and then will be added to the dataset.
            for f = 1 : length(Fracture)
                for If = 1 : Fracture(f).N_Total
                    for g = 1 : length(Fracture)
                        if f ~= g && ~isempty(Fracture(f).CI_frac_pEDFM{g})
                            if ~isempty(Fracture(f).CI_frac_pEDFM{g}{If})
                                Ig_List = Fracture(f).CI_frac_pEDFM{g}{If}(:,1);
                                
                                for n = 1 : length(Ig_List)
                                    bothExist = false;
                                    Ig = Ig_List(n);
                                    CI_frac1 = Fracture(f).CI_frac_pEDFM{g}{If}(n,2);
                                    if ~isempty( Fracture(g).CI_frac_pEDFM{f}{Ig} )
                                        ind_List = Fracture(g).CI_frac_pEDFM{f}{Ig}(:,1);
                                        if ismember(If,ind_List)
                                            ind = find(ind_List == If);
                                            CI_frac2 = Fracture(g).CI_frac_pEDFM{f}{Ig}(ind,2);
                                            bothExist = true;
                                        else
                                            CI_frac2 = CI_frac1;
                                        end
                                    else
                                        CI_frac2 = CI_frac1;
                                    end
                                    CI_Harmonic = 2 * CI_frac1 * CI_frac2 / (CI_frac1 + CI_frac2);
                                    Fracture(f).CI_frac_pEDFM{g}{If}(n,2) = CI_Harmonic;
                                    if ~bothExist
                                        Fracture(g).CI_frac_pEDFM{f}{Ig} = [ Fracture(g).CI_frac_pEDFM{f}{Ig} ; [If , CI_Harmonic] ];
                                        Fracture(g).CI_frac_pEDFM{f}{Ig} = sortrows(Fracture(g).CI_frac_pEDFM{f}{Ig});
                                        Fracture(g).NumOf_fracConn_pEDFM(Ig) = Fracture(g).NumOf_fracConn_pEDFM(Ig) + 1;
                                    else
                                        Fracture(g).CI_frac_pEDFM{f}{Ig}(ind,2) = CI_Harmonic;
                                    end
                                end
                                
                            end
                        end
                    end
                end
            end
            
            % Improving pEDFM alpha corrections
            obj.ReservoirGrid.alpha_Tx( obj.ReservoirGrid.alpha_Tx > 0.9 ) = 1;
            obj.ReservoirGrid.alpha_Ty( obj.ReservoirGrid.alpha_Ty > 0.9 ) = 1;
            obj.ReservoirGrid.alpha_Tz( obj.ReservoirGrid.alpha_Tz > 0.9 ) = 1;
            for f = 1 : length(Fracture)
                ind = find( Fracture(f).alpha_Tx > 0 );
                if length(ind) == 1
                    Fracture(f).alpha_Tx(ind) = 1;
                end
                Fracture(f).alpha_Tx( Fracture(f).alpha_Tx > 0.9 ) = 1;
                ind = find( Fracture(f).alpha_Ty > 0 );
                if length(ind) == 1
                    Fracture(f).alpha_Ty(ind) = 1;
                end
                Fracture(f).alpha_Ty( Fracture(f).alpha_Ty > 0.9 ) = 1;
            end
            
            fprintf('\n---------------------------------------------------------\n');
        end
        function PlotFractures(obj)
            %% Plotting fracture plates
            Fracture = obj.FracturesGrid.Fracture;
            fprintf('Plotting fractures ...\n');
            figure();
            for f = 1 : length(Fracture)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Plotting the boundaries of each fracture plate
                plot3(Fracture(f).Points(1,[1:4,1]),Fracture(f).Points(2,[1:4,1]),Fracture(f).Points(3,[1:4,1]),'r','LineWidth',3);
                xlim([0,obj.ReservoirGrid.LX]); ylim([0,obj.ReservoirGrid.LY]); zlim([0,obj.ReservoirGrid.LZ]);
                xlabel('X [cm]'); ylabel('Y [cm]'); zlabel('Z [cm]');
                grid ON;
                
                % Plotting the normal vector of each fracture plate
                hold on;
                quiver3(Fracture(f).PointM(1),Fracture(f).PointM(2),Fracture(f).PointM(3),Fracture(f).n_vec_plot(1),Fracture(f).n_vec_plot(2),Fracture(f).n_vec_plot(3),'b','LineWidth',2);
                
                % Plotting cell edges of each fracture plate
                for i = 1 : Fracture(f).N_Length_AB+1
                    plot3(Fracture(f).GridCoords(i,:,1),Fracture(f).GridCoords(i,:,2),Fracture(f).GridCoords(i,:,3),'r','LineWidth',1);
                end
                for j = 1 : Fracture(f).N_Width_AD+1
                    plot3(Fracture(f).GridCoords(:,j,1),Fracture(f).GridCoords(:,j,2),Fracture(f).GridCoords(:,j,3),'r','LineWidth',1);
                end
                
                % Writing the letters of A,B,C,D to the corners
                text( Fracture(f).PointA(1) - 0.05*Fracture(f).AD_vec(1) - 0.05*Fracture(f).AB_vec(1) , ...
                      Fracture(f).PointA(2) - 0.05*Fracture(f).AD_vec(2) - 0.05*Fracture(f).AB_vec(2) , ...
                      Fracture(f).PointA(3) - 0.05*Fracture(f).AD_vec(3) - 0.05*Fracture(f).AB_vec(3) , ...
                      'A' , 'Color' , 'black' , 'FontSize' , 14);
                text( Fracture(f).PointB(1) - 0.05*Fracture(f).AD_vec(1) + 0.05*Fracture(f).AB_vec(1) , ...
                      Fracture(f).PointB(2) - 0.05*Fracture(f).AD_vec(2) + 0.05*Fracture(f).AB_vec(2) , ...
                      Fracture(f).PointB(3) - 0.05*Fracture(f).AD_vec(3) + 0.05*Fracture(f).AB_vec(3) , ...
                      'B' , 'Color' , 'black' , 'FontSize' , 14);
                text( Fracture(f).PointC(1) + 0.05*Fracture(f).AD_vec(1) + 0.05*Fracture(f).AB_vec(1) , ...
                      Fracture(f).PointC(2) + 0.05*Fracture(f).AD_vec(2) + 0.05*Fracture(f).AB_vec(2) , ...
                      Fracture(f).PointC(3) + 0.05*Fracture(f).AD_vec(3) + 0.05*Fracture(f).AB_vec(3) , ...
                      'C' , 'Color' , 'black' , 'FontSize' , 14);
                text( Fracture(f).PointD(1) + 0.05*Fracture(f).AD_vec(1) - 0.05*Fracture(f).AB_vec(1) , ...
                      Fracture(f).PointD(2) + 0.05*Fracture(f).AD_vec(2) - 0.05*Fracture(f).AB_vec(2) , ...
                      Fracture(f).PointD(3) + 0.05*Fracture(f).AD_vec(3) - 0.05*Fracture(f).AB_vec(3) , ...
                      'D' , 'Color' , 'black' , 'FontSize' , 14);
                text( Fracture(f).PointM(1), Fracture(f).PointM(2), Fracture(f).PointM(3), ...
                      num2str(f) , 'Color' , 'magenta' , 'FontSize' , 20, 'FontWeight','bold');
                
                % Plotting the intersection line of each two fracture plates
                for g = f+1 : length(Fracture)
                    if ~isempty(Fracture(f).intersectCoord_fracturePlate{g})
                        plot3(Fracture(f).intersectCoord_fracturePlate{g}(1,:),...
                              Fracture(f).intersectCoord_fracturePlate{g}(2,:),...
                              Fracture(f).intersectCoord_fracturePlate{g}(3,:),...
                              'Color',[0,0.5,0],'LineWidth',5);
                    end
                end
            end
            
            rotate3d on;
            view([0 90]);
            fprintf('---------------------------------------------------------\n');
        end
        function [intersectCoordFinal, intersectCoordTemp , areCoplanar] = Obtain_PlaneSegment_Cube_Intersection(obj, i,j,k, Plane_fracCell)
            % Assigning reservoir properties
            % Coordinates of reservoir grid nodes in each direction
            Xim = obj.ReservoirGrid.Xim;
            Yim = obj.ReservoirGrid.Yim;
            Zim = obj.ReservoirGrid.Zim;
            Epsilon = obj.ReservoirGrid.Epsilon;
            
            % Initializing some variables
            doNotContinue       = 0;
            areCoplanar         = 0;
            intersectCoordTemp  = []; 
            intersectCoordFinal = [];
            
            % Matrix Cube Face #X1
            if doNotContinue == 0
                PointA = [ Xim(i  ) ; Yim(j  ) ; Zim(k  ) ];
                PointB = [ Xim(i  ) ; Yim(j+1) ; Zim(k  ) ];
                PointC = [ Xim(i  ) ; Yim(j+1) ; Zim(k+1) ];
                PointD = [ Xim(i  ) ; Yim(j  ) ; Zim(k+1) ];
                Plane_matFace_X1 = planeSegment_FracGen(PointA,PointB,PointC,PointD);
                
                [Geostatus_X1, intersectPoints_X1] = Plane_fracCell.Obtain_PlaneSegment_PlaneSegment_Intersection( Plane_matFace_X1 , Epsilon );
                if ( Geostatus_X1.areCoplanar == 1 )
                    doNotContinue      = 1;
                    areCoplanar        = 1;
                    % intersectPoints_X1 = [ intersectPoints_X1 , Plane_matFace_X1.PointA , Plane_matFace_X1.PointB , Plane_matFace_X1.PointC , Plane_matFace_X1.PointD ];
                end
                intersectCoordTemp = [ intersectCoordTemp , intersectPoints_X1 ];
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Matrix Cube Face #X2
            if doNotContinue == 0
                PointA = [ Xim(i+1) ; Yim(j  ) ; Zim(k  ) ];
                PointB = [ Xim(i+1) ; Yim(j+1) ; Zim(k  ) ];
                PointC = [ Xim(i+1) ; Yim(j+1) ; Zim(k+1) ];
                PointD = [ Xim(i+1) ; Yim(j  ) ; Zim(k+1) ];
                Plane_matFace_X2 = planeSegment_FracGen(PointA,PointB,PointC,PointD);
                
                [Geostatus_X2, intersectPoints_X2] = Plane_fracCell.Obtain_PlaneSegment_PlaneSegment_Intersection( Plane_matFace_X2 , Epsilon );
                if ( Geostatus_X2.areCoplanar == 1 )
                    doNotContinue      = 1;
                    areCoplanar        = 1;
                    % intersectPoints_X2 = [ intersectPoints_X2 , Plane_matFace_X2.PointA , Plane_matFace_X2.PointB , Plane_matFace_X2.PointC , Plane_matFace_X2.PointD ];
                end
                intersectCoordTemp = [ intersectCoordTemp , intersectPoints_X2 ];
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Matrix Cube Face #Y1
            if doNotContinue == 0
                PointA = [ Xim(i  ) ; Yim(j  ) ; Zim(k  ) ];
                PointB = [ Xim(i+1) ; Yim(j  ) ; Zim(k  ) ];
                PointC = [ Xim(i+1) ; Yim(j  ) ; Zim(k+1) ];
                PointD = [ Xim(i  ) ; Yim(j  ) ; Zim(k+1) ];
                Plane_matFace_Y1 = planeSegment_FracGen(PointA,PointB,PointC,PointD);
                
                [Geostatus_Y1, intersectPoints_Y1] = Plane_fracCell.Obtain_PlaneSegment_PlaneSegment_Intersection( Plane_matFace_Y1 , Epsilon );
                if ( Geostatus_Y1.areCoplanar == 1 )
                    doNotContinue      = 1;
                    areCoplanar        = 1;
                    % intersectPoints_Y1 = [ intersectPoints_Y1 , Plane_matFace_Y1.PointA , Plane_matFace_Y1.PointB , Plane_matFace_Y1.PointC , Plane_matFace_Y1.PointD ];
                end
                intersectCoordTemp = [ intersectCoordTemp , intersectPoints_Y1 ];
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Matrix Cube Face #Y2
            if doNotContinue == 0
                PointA = [ Xim(i  ) ; Yim(j+1) ; Zim(k  ) ];
                PointB = [ Xim(i+1) ; Yim(j+1) ; Zim(k  ) ];
                PointC = [ Xim(i+1) ; Yim(j+1) ; Zim(k+1) ];
                PointD = [ Xim(i  ) ; Yim(j+1) ; Zim(k+1) ];
                Plane_matFace_Y2 = planeSegment_FracGen(PointA,PointB,PointC,PointD);
                
                [Geostatus_Y2, intersectPoints_Y2] = Plane_fracCell.Obtain_PlaneSegment_PlaneSegment_Intersection( Plane_matFace_Y2 , Epsilon );
                if ( Geostatus_Y2.areCoplanar == 1 )
                    doNotContinue      = 1;
                    areCoplanar        = 1;
                    % intersectPoints_Y2 = [ intersectPoints_Y2 , Plane_matFace_Y2.PointA , Plane_matFace_Y2.PointB , Plane_matFace_Y2.PointC , Plane_matFace_Y2.PointD ];
                end
                intersectCoordTemp = [ intersectCoordTemp , intersectPoints_Y2 ];
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Matrix Cube Face #Z1
            if doNotContinue == 0
                PointA = [ Xim(i  ) ; Yim(j  ) ; Zim(k  ) ];
                PointB = [ Xim(i+1) ; Yim(j  ) ; Zim(k  ) ];
                PointC = [ Xim(i+1) ; Yim(j+1) ; Zim(k  ) ];
                PointD = [ Xim(i  ) ; Yim(j+1) ; Zim(k  ) ];
                Plane_matFace_Z1 = planeSegment_FracGen(PointA,PointB,PointC,PointD);
                
                [Geostatus_Z1, intersectPoints_Z1] = Plane_fracCell.Obtain_PlaneSegment_PlaneSegment_Intersection( Plane_matFace_Z1 , Epsilon );
                if ( Geostatus_Z1.areCoplanar == 1 )
                    doNotContinue      = 1;
                    areCoplanar        = 1;
                    % intersectPoints_Z1 = [ intersectPoints_Z1 , Plane_matFace_Z1.PointA , Plane_matFace_Z1.PointB , Plane_matFace_Z1.PointC , Plane_matFace_Z1.PointD ];
                end
                intersectCoordTemp = [ intersectCoordTemp , intersectPoints_Z1 ];
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Matrix Cube Face #Z2
            if doNotContinue == 0
                PointA = [ Xim(i  ) ; Yim(j  ) ; Zim(k+1) ];
                PointB = [ Xim(i+1) ; Yim(j  ) ; Zim(k+1) ];
                PointC = [ Xim(i+1) ; Yim(j+1) ; Zim(k+1) ];
                PointD = [ Xim(i  ) ; Yim(j+1) ; Zim(k+1) ];
                Plane_matFace_Z2 = planeSegment_FracGen(PointA,PointB,PointC,PointD);
                
                [Geostatus_Z2, intersectPoints_Z2] = Plane_fracCell.Obtain_PlaneSegment_PlaneSegment_Intersection( Plane_matFace_Z2 , Epsilon );
                if ( Geostatus_Z2.areCoplanar == 1 )
                    doNotContinue = 1;
                    areCoplanar   = 1;
                    % intersectPoints_Z2 = [ intersectPoints_Z2 , Plane_matFace_Z2.PointA , Plane_matFace_Z2.PointB , Plane_matFace_Z2.PointC , Plane_matFace_Z2.PointD ];
                end
                intersectCoordTemp = [ intersectCoordTemp , intersectPoints_Z2 ];
            end
            
            intersectCoordTemp = [ intersectCoordTemp , Plane_fracCell.PointA , Plane_fracCell.PointB , Plane_fracCell.PointC , Plane_fracCell.PointD ];
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Processing the intersection Points
            for nr = 1 : size( intersectCoordTemp , 2 )
                
                % Removing the Points that are not inside the matrix cell
                if ( ( Xim(i) - intersectCoordTemp(1,nr) ) > Epsilon ) || ( ( intersectCoordTemp(1,nr) - Xim(i+1) ) > Epsilon )
                    intersectCoordTemp(:,nr) = [ NaN ; NaN ; NaN ];
                    continue;
                end
                if ( ( Yim(j) - intersectCoordTemp(2,nr) ) > Epsilon ) || ( ( intersectCoordTemp(2,nr) - Yim(j+1) ) > Epsilon )
                    intersectCoordTemp(:,nr) = [ NaN ; NaN ; NaN ];
                    continue;
                end
                if ( ( Zim(k) - intersectCoordTemp(3,nr) ) > Epsilon ) || ( ( intersectCoordTemp(3,nr) - Zim(k+1) ) > Epsilon )
                    intersectCoordTemp(:,nr) = [ NaN ; NaN ; NaN ];
                    continue;
                end
                
                % Removing the Points that are not inside the fracture cell
                isInside = Plane_fracCell.Is_Point_Inside_PlaneSegment(intersectCoordTemp(:,nr), Epsilon);
                if isInside == 0
                    intersectCoordTemp(:,nr) = [ NaN ; NaN ; NaN ];
                    continue;
                end
                
                % Removing the Points that are repeated
                subtracted = ( intersectCoordTemp - intersectCoordTemp(:,nr) );
                subtracted = sqrt( subtracted(1,:).^2 + subtracted(2,:).^2 + subtracted(3,:).^2 );
                if ~isempty( find( subtracted([1:nr-1 , nr+1:end]) <= Epsilon , 1 ) )
                    intersectCoordTemp(:,nr) = [ NaN ; NaN ; NaN ];
                    continue;
                end
                
                % Writing the intersection Points into a new variable
                if ~isnan(intersectCoordTemp(1,nr))
                    intersectCoordFinal = [ intersectCoordFinal , intersectCoordTemp(:,nr) ];
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Ordering the intersection Points (either clockwise or counter clockwise)
            intersectCoordTemp = intersectCoordFinal;
            if size(intersectCoordTemp,2) > 2
                distanceToZero = sqrt( sum( intersectCoordTemp.^2 , 1 ) );
                
                % Setting the 1st Point
                [ ~ , index ] = min( distanceToZero );
                intersectCoordFinal      = zeros ( size(intersectCoordTemp,1) , size(intersectCoordTemp,2) );
                intersectCoordFinal(:,1) = intersectCoordTemp(:,index);
                intersectCoordTemp(:,index) = [];
                
                % Setting the 2nd Point
                distanceToRef = sqrt( sum ( ( intersectCoordTemp - intersectCoordFinal(:,1) ).^2 , 1 ) );
                [ ~ , index ] = min( distanceToRef );
                intersectCoordFinal(:,2) = intersectCoordTemp(:,index);
                intersectCoordTemp(:,index) = [];
                
                % Setting the rest of the Points based on the maximum angle between two close segments
                for m = 2 : length( distanceToZero ) - 2
                    arcCos_Theta = ones( size( intersectCoordTemp , 2 ) , 1 );
                    for n = 1 : size( intersectCoordTemp , 2 )
                        arcCos_Theta(n) = dot( intersectCoordFinal(:,m) - intersectCoordFinal(:,m-1) , ...
                            intersectCoordFinal(:,m) - intersectCoordTemp (:,n  ) ) ...
                            / norm( intersectCoordFinal(:,m) - intersectCoordFinal(:,m-1) ) ...
                            / norm( intersectCoordFinal(:,m) - intersectCoordTemp (:,n  ) );
                    end
                    [ ~ , index ] = min( arcCos_Theta );
                    intersectCoordFinal(:,m+1)  = intersectCoordTemp(:,index);
                    intersectCoordTemp(:,index) = [];
                end
                intersectCoordFinal(:,end) = intersectCoordTemp;
            end
        end
        function [ Af , AvgDist , Collinearity] = Line_Plane_Connectivity_3D( obj, Plane , Line )
            Epsilon = obj.ReservoirGrid.Epsilon;
            Collinearity = 0;
            
            %% Calculating Af of the line segment inside the plane segment
            Af = 2 * norm ( Line.AB_vec );
            AB_Inf = lineInf_FracGen(Plane.AB_vec);  AB_Inf.AddPoint0(Plane.PointA);  AB_Inf.AddEquation(Plane.PointA);
            BC_Inf = lineInf_FracGen(Plane.BC_vec);  BC_Inf.AddPoint0(Plane.PointB);  BC_Inf.AddEquation(Plane.PointA);
            CD_Inf = lineInf_FracGen(Plane.CD_vec);  CD_Inf.AddPoint0(Plane.PointC);  CD_Inf.AddEquation(Plane.PointA);
            DA_Inf = lineInf_FracGen(Plane.DA_vec);  DA_Inf.AddPoint0(Plane.PointD);  DA_Inf.AddEquation(Plane.PointA);
            
            [Geostatus_AB, IntersectPoint_AB] = Line.Obtain_LineInf_LineInf_Intersection( AB_Inf, Epsilon );
            [Geostatus_BC, IntersectPoint_BC] = Line.Obtain_LineInf_LineInf_Intersection( BC_Inf, Epsilon );
            [Geostatus_CD, IntersectPoint_CD] = Line.Obtain_LineInf_LineInf_Intersection( CD_Inf, Epsilon );
            [Geostatus_DA, IntersectPoint_DA] = Line.Obtain_LineInf_LineInf_Intersection( DA_Inf, Epsilon );
            
            if ( Geostatus_AB.areCollinear == 1 ) || ( Geostatus_BC.areCollinear == 1 ) || ...
               ( Geostatus_CD.areCollinear == 1 ) || ( Geostatus_DA.areCollinear == 1 )
                Collinearity = 1;
                Af = Af /2;
            end
            
            %% Obtaining aveDist of the plane segment from the line segment
            AvgDist  = Plane.Obtain_Average_Dsitance_LineInf_From_PlaneSegment(Line, 6);
            
%             % Parallel situation
%             if ( Geostatus_AB.areParallel == 1 ) || ( Geostatus_CD.areParallel == 1 )
%                 % The line segment is parallel to AB and CD
%                 Dist1 = norm( cross( (Line.PointA - Plane.PointA) , Line.AB_vec ) ) / norm(Line.AB_vec);
%                 Dist2 = norm( cross( (Line.PointA - Plane.PointC) , Line.AB_vec ) ) / norm(Line.AB_vec);
%                 if abs(Dist1+Dist2-norm(Plane.DA_vec)) > Epsilon , error('The summation of distances of the parallel line from AB and CD is not correct!');  end
%                 aveDist = ( Dist1^2 + Dist2^2 ) / ( 2*Dist1 + 2*Dist2 );
%                 return;
%             end
%             if ( Geostatus_BC.areParallel == 1 ) || ( Geostatus_DA.areParallel == 1 )
%                 % The line segment is parallel to BC and DA
%                 Dist1 = norm( cross( (Line.PointA - Plane.PointB) , Line.AB_vec ) ) / norm(Line.AB_vec);
%                 Dist2 = norm( cross( (Line.PointA - Plane.PointD) , Line.AB_vec ) ) / norm(Line.AB_vec);
%                 if abs(Dist1+Dist2-norm(Plane.AB_vec)) > Epsilon , error('The summation of distances of the parallel line from BC and DA is not correct!');  end
%                 aveDist = ( Dist1^2 + Dist2^2 ) / ( 2*Dist1 + 2*Dist2 );
%                 return;
%             end
%             
%             % Non-parallel situation
%             isInside_AB = Plane.Is_Point_Inside_PlaneSegment(IntersectPoint_AB, Epsilon);
%             isInside_BC = Plane.Is_Point_Inside_PlaneSegment(IntersectPoint_BC, Epsilon);
%             isInside_CD = Plane.Is_Point_Inside_PlaneSegment(IntersectPoint_CD, Epsilon);
%             isInside_DA = Plane.Is_Point_Inside_PlaneSegment(IntersectPoint_DA, Epsilon);
%             
%             if sum( [ isInside_AB , isInside_BC , isInside_CD , isInside_DA ] ) == 4
%                 % The line segment is on diagonal of plane segment
%                 aveDist = norm(Plane.AB_vec) * norm(Plane.BC_vec) / ( 3 * sqrt( norm(Plane.AB_vec)^2 + norm(Plane.BC_vec)^2 ) );
%                 return;
%             end
%             
%             if ( sum( [ isInside_AB , isInside_CD ] ) == 2 )
%                 % The extension of the line segment intersects with two parallel sides of the plane segment (AB and CD)
%                 if norm( Plane.PointA - IntersectPoint_DA ) < norm( Plane.PointD - IntersectPoint_DA )
%                     Area1 = Triangle_Area_3D( IntersectPoint_DA , Plane.PointD , IntersectPoint_CD );  % Large Triangle
%                     Area2 = Triangle_Area_3D( IntersectPoint_DA , Plane.PointA , IntersectPoint_AB );  % Small Triangle
%                     Area3 = Triangle_Area_3D( IntersectPoint_BC , Plane.PointB , IntersectPoint_AB );  % Large Triangle
%                     Area4 = Triangle_Area_3D( IntersectPoint_BC , Plane.PointC , IntersectPoint_CD );  % Small Triangle
%                     Length1 = norm( IntersectPoint_CD - Plane.PointD );   Width1 = norm( Plane.PointD - IntersectPoint_CD );  % Large Triangle
%                     Length2 = norm( IntersectPoint_AB - Plane.PointA );   Width2 = norm( Plane.PointA - IntersectPoint_AB );  % Small Triangle
%                     Length3 = norm( IntersectPoint_AB - Plane.PointB );   Width3 = norm( Plane.PointB - IntersectPoint_AB );  % Large Triangle
%                     Length4 = norm( IntersectPoint_CD - Plane.PointC );   Width4 = norm( Plane.PointC - IntersectPoint_CD );  % Small Triangle
%                 else
%                     Area1 = Triangle_Area_3D( IntersectPoint_DA , Plane.PointA , IntersectPoint_AB );  % Large Triangle
%                     Area2 = Triangle_Area_3D( IntersectPoint_DA , Plane.PointD , IntersectPoint_CD );  % Small Triangle
%                     Area3 = Triangle_Area_3D( IntersectPoint_BC , Plane.PointC , IntersectPoint_CD );  % Large Triangle
%                     Area4 = Triangle_Area_3D( IntersectPoint_BC , Plane.PointB , IntersectPoint_AB );  % Small Triangle
%                     Length1 = norm( IntersectPoint_AB - Plane.PointA );   Width1 = norm( Plane.PointA - IntersectPoint_AB );  % Large Triangle
%                     Length2 = norm( IntersectPoint_CD - Plane.PointD );   Width2 = norm( Plane.PointD - IntersectPoint_CD );  % Small Triangle
%                     Length3 = norm( IntersectPoint_CD - Plane.PointC );   Width3 = norm( Plane.PointC - IntersectPoint_CD );  % Large Triangle
%                     Length4 = norm( IntersectPoint_AB - Plane.PointB );   Width4 = norm( Plane.PointB - IntersectPoint_AB );  % Small Triangle
%                 end
%                 Dist1 = Length1 * Width1 / ( 3 * sqrt( Length1^2 + Width1^2 ) );  if Area1==0,  Dist1=0;  end
%                 Dist2 = Length2 * Width2 / ( 3 * sqrt( Length2^2 + Width2^2 ) );  if Area2==0,  Dist2=0;  end
%                 Dist3 = Length3 * Width3 / ( 3 * sqrt( Length3^2 + Width3^2 ) );  if Area3==0,  Dist3=0;  end
%                 Dist4 = Length4 * Width4 / ( 3 * sqrt( Length4^2 + Width4^2 ) );  if Area4==0,  Dist4=0;  end
%                 aveDist = ( Area1*Dist1 - Area2*Dist2 + Area3*Dist3 - Area4*Dist4 ) / ( Area1 - Area2 + Area3 - Area4 );
%                 return;
%             end
%             
%             if ( sum( [ isInside_BC , isInside_DA ] ) == 2 )
%                 % The extension of the line segment intersects with two parallel sides of the plane segment (BC and DA)
%                 if norm( Plane.PointA - IntersectPoint_AB ) < norm( Plane.PointB - IntersectPoint_AB )
%                     Area1 = Triangle_Area_3D( IntersectPoint_AB , Plane.PointB , IntersectPoint_BC );  % Large Triangle
%                     Area2 = Triangle_Area_3D( IntersectPoint_AB , Plane.PointA , IntersectPoint_DA );  % Small Triangle
%                     Area3 = Triangle_Area_3D( IntersectPoint_CD , Plane.PointD , IntersectPoint_DA );  % Large Triangle
%                     Area4 = Triangle_Area_3D( IntersectPoint_CD , Plane.PointC , IntersectPoint_BC );  % Small Triangle
%                     Length1 = norm( IntersectPoint_AB - Plane.PointB );   Width1 = norm( Plane.PointB - IntersectPoint_BC );  % Large Triangle
%                     Length2 = norm( IntersectPoint_AB - Plane.PointA );   Width2 = norm( Plane.PointA - IntersectPoint_DA );  % Small Triangle
%                     Length3 = norm( IntersectPoint_CD - Plane.PointD );   Width3 = norm( Plane.PointD - IntersectPoint_DA );  % Large Triangle
%                     Length4 = norm( IntersectPoint_CD - Plane.PointC );   Width4 = norm( Plane.PointC - IntersectPoint_BC );  % Small Triangle
%                 else
%                     Area1 = Triangle_Area_3D( IntersectPoint_AB , Plane.PointA , IntersectPoint_DA );  % Large Triangle
%                     Area2 = Triangle_Area_3D( IntersectPoint_AB , Plane.PointB , IntersectPoint_BC );  % Small Triangle
%                     Area3 = Triangle_Area_3D( IntersectPoint_CD , Plane.PointC , IntersectPoint_BC );  % Large Triangle
%                     Area4 = Triangle_Area_3D( IntersectPoint_CD , Plane.PointD , IntersectPoint_DA );  % Small Triangle
%                     Length1 = norm( IntersectPoint_AB - Plane.PointA );   Width1 = norm( Plane.PointA - IntersectPoint_DA );  % Large Triangle
%                     Length2 = norm( IntersectPoint_AB - Plane.PointB );   Width2 = norm( Plane.PointB - IntersectPoint_BC );  % Small Triangle
%                     Length3 = norm( IntersectPoint_CD - Plane.PointC );   Width3 = norm( Plane.PointC - IntersectPoint_BC );  % Large Triangle
%                     Length4 = norm( IntersectPoint_CD - Plane.PointD );   Width4 = norm( Plane.PointD - IntersectPoint_DA );  % Small Triangle
%                 end
%                 Dist1 = Length1 * Width1 / ( 3 * sqrt( Length1^2 + Width1^2 ) );  if Area1==0,  Dist1=0;  end
%                 Dist2 = Length2 * Width2 / ( 3 * sqrt( Length2^2 + Width2^2 ) );  if Area2==0,  Dist2=0;  end
%                 Dist3 = Length3 * Width3 / ( 3 * sqrt( Length3^2 + Width3^2 ) );  if Area3==0,  Dist3=0;  end
%                 Dist4 = Length4 * Width4 / ( 3 * sqrt( Length4^2 + Width4^2 ) );  if Area4==0,  Dist4=0;  end
%                 aveDist = ( Area1*Dist1 - Area2*Dist2 + Area3*Dist3 - Area4*Dist4 ) / ( Area1 - Area2 + Area3 - Area4 );
%                 return;
%             end
%             
%             if ( sum( [ isInside_AB , isInside_BC ] ) == 2 )
%                 % The extension of the line segment intersects with two cornering sides of the plane segment (AB and BC)
%                 Area1 = Triangle_Area_3D( IntersectPoint_CD , Plane.PointD , IntersectPoint_DA );  % Largest Triangle
%                 Area2 = Triangle_Area_3D( IntersectPoint_AB , Plane.PointB , IntersectPoint_BC );  % Small Triangle Inside
%                 Area3 = Triangle_Area_3D( IntersectPoint_DA , Plane.PointA , IntersectPoint_AB );  % Outer Triangle 1
%                 Area4 = Triangle_Area_3D( IntersectPoint_BC , Plane.PointC , IntersectPoint_CD );  % Outer Triangle 2
%                 Length1 = norm( IntersectPoint_CD - Plane.PointD );   Width1 = norm( Plane.PointD - IntersectPoint_DA );  % Large Triangle
%                 Length2 = norm( IntersectPoint_AB - Plane.PointB );   Width2 = norm( Plane.PointB - IntersectPoint_BC );  % Small Triangle
%                 Length3 = norm( IntersectPoint_DA - Plane.PointA );   Width3 = norm( Plane.PointA - IntersectPoint_AB );  % Large Triangle
%                 Length4 = norm( IntersectPoint_BC - Plane.PointC );   Width4 = norm( Plane.PointC - IntersectPoint_CD );  % Small Triangle
%                 
%             elseif ( sum( [ isInside_BC , isInside_CD ] ) == 2 )
%                 % The extension of the line segment intersects with two cornering sides of the plane segment (BC and CD)
%                 Area1 = Triangle_Area_3D( IntersectPoint_DA , Plane.PointA , IntersectPoint_AB );  % Largest Triangle
%                 Area2 = Triangle_Area_3D( IntersectPoint_BC , Plane.PointC , IntersectPoint_CD );  % Small Triangle Inside
%                 Area3 = Triangle_Area_3D( IntersectPoint_AB , Plane.PointB , IntersectPoint_BC );  % Outer Triangle 1
%                 Area4 = Triangle_Area_3D( IntersectPoint_CD , Plane.PointD , IntersectPoint_DA );  % Outer Triangle 2
%                 Length1 = norm( IntersectPoint_DA - Plane.PointA );   Width1 = norm( Plane.PointA - IntersectPoint_AB );  % Large Triangle
%                 Length2 = norm( IntersectPoint_BC - Plane.PointC );   Width2 = norm( Plane.PointC - IntersectPoint_CD );  % Small Triangle
%                 Length3 = norm( IntersectPoint_AB - Plane.PointB );   Width3 = norm( Plane.PointB - IntersectPoint_BC );  % Large Triangle
%                 Length4 = norm( IntersectPoint_CD - Plane.PointD );   Width4 = norm( Plane.PointD - IntersectPoint_DA );  % Small Triangle
%                 
%             elseif ( sum( [ isInside_CD , isInside_DA ] ) == 2 )
%                 % The extension of the line segment intersects with two cornering sides of the plane segment (CD and DA)
%                 Area1 = Triangle_Area_3D( IntersectPoint_AB , Plane.PointB , IntersectPoint_BC );  % Largest Triangle
%                 Area2 = Triangle_Area_3D( IntersectPoint_CD , Plane.PointD , IntersectPoint_DA );  % Small Triangle Inside
%                 Area3 = Triangle_Area_3D( IntersectPoint_DA , Plane.PointA , IntersectPoint_AB );  % Outer Triangle 1
%                 Area4 = Triangle_Area_3D( IntersectPoint_BC , Plane.PointC , IntersectPoint_CD );  % Outer Triangle 2
%                 Length1 = norm( IntersectPoint_AB - Plane.PointB );   Width1 = norm( Plane.PointB - IntersectPoint_BC );  % Large Triangle
%                 Length2 = norm( IntersectPoint_CD - Plane.PointD );   Width2 = norm( Plane.PointD - IntersectPoint_DA );  % Small Triangle
%                 Length3 = norm( IntersectPoint_DA - Plane.PointA );   Width3 = norm( Plane.PointA - IntersectPoint_AB );  % Large Triangle
%                 Length4 = norm( IntersectPoint_BC - Plane.PointC );   Width4 = norm( Plane.PointC - IntersectPoint_CD );  % Small Triangle
%                 
%             elseif ( sum( [ isInside_DA , isInside_AB ] ) == 2 )
%                 % The extension of the line segment intersects with two cornering sides of the plane segment (DA and AB)
%                 Area1 = Triangle_Area_3D( IntersectPoint_BC , Plane.PointC , IntersectPoint_CD );  % Largest Triangle
%                 Area2 = Triangle_Area_3D( IntersectPoint_DA , Plane.PointA , IntersectPoint_AB );  % Small Triangle Inside
%                 Area3 = Triangle_Area_3D( IntersectPoint_AB , Plane.PointB , IntersectPoint_BC );  % Outer Triangle 1
%                 Area4 = Triangle_Area_3D( IntersectPoint_CD , Plane.PointD , IntersectPoint_DA );  % Outer Triangle 2
%                 Length1 = norm( IntersectPoint_BC - Plane.PointC );   Width1 = norm( Plane.PointC - IntersectPoint_CD );  % Large Triangle
%                 Length2 = norm( IntersectPoint_DA - Plane.PointA );   Width2 = norm( Plane.PointA - IntersectPoint_AB );  % Small Triangle
%                 Length3 = norm( IntersectPoint_AB - Plane.PointB );   Width3 = norm( Plane.PointB - IntersectPoint_BC );  % Large Triangle
%                 Length4 = norm( IntersectPoint_CD - Plane.PointD );   Width4 = norm( Plane.PointD - IntersectPoint_DA );  % Small Triangle
%             end
%             
%             Dist1 = Length1 * Width1 / ( 3 * sqrt( Length1^2 + Width1^2 ) );  if Area1==0,  Dist1=0;  end
%             Dist2 = Length2 * Width2 / ( 3 * sqrt( Length2^2 + Width2^2 ) );  if Area2==0,  Dist2=0;  end
%             Dist3 = Length3 * Width3 / ( 3 * sqrt( Length3^2 + Width3^2 ) );  if Area3==0,  Dist3=0;  end
%             Dist4 = Length4 * Width4 / ( 3 * sqrt( Length4^2 + Width4^2 ) );  if Area4==0,  Dist4=0;  end
%             aveDist = ( Area1*Dist1 + Area2*Dist2 - Area3*Dist3 - Area4*Dist4 ) / ( Area1 + Area2 - Area3 - Area4 );
        end
    end
end