%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FracGen2D  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Mousa HosseiniMehr, MSc Petroleum Engineering, CEG Faculty, TU Delft
% Project: 3D EDFM Package for F-ADM, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build Date : 2016-11-29
% Modified on: 2017-03-17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Frac = FracGen2D(obj)

LX = obj.Simulation.Reservoir.LX;
LY = obj.Simulation.Reservoir.LY;
LZ = obj.Simulation.Reservoir.LZ;
NX = obj.Simulation.Reservoir.NX;
NY = obj.Simulation.Reservoir.NY;
NZ = obj.Simulation.Reservoir.NZ;
Xcm = obj.Simulation.Reservoir.Xcm;
Ycm = obj.Simulation.Reservoir.Ycm;
Zcm = obj.Simulation.Reservoir.Zcm;
Xim = obj.Simulation.Reservoir.Xim;
Yim = obj.Simulation.Reservoir.Yim;
Zim = obj.Simulation.Reservoir.Zim;

Frac_Input = obj.Builder.FracInput;

almostZero = obj.Simulation.Reservoir.almostZero;

Frac = fractures_mousa();
%% Fractures Construction with Geometry Input
f = 0;
for F = 1 : size(Frac_Input,2)
    if Frac_Input(F).isActive == 1
        f = f + 1;
    
        if isnan( Frac_Input(F).CenterCoord(1) )
            % Fractures construction with 2 Points input type

            % Setting attributes to the corners (A is the corner with 90 deg angle)
            Frac(f).PointA = Frac_Input(F).CornerCoords(:,1);  Frac(f).PointA(3)=0;
            Frac(f).PointB = Frac_Input(F).CornerCoords(:,2);  Frac(f).PointB(3)=0;

            % The centeral Point of the line 
            Frac(f).PointM = ( Frac(f).PointA + Frac(f).PointB ) /2;

            % The coordinates of all Points in one variable (matrix) for Plotting
            Frac(f).Points = [ Frac(f).PointA , Frac(f).PointB , Frac(f).PointM ];

        else
            % Fracture construction with central Point and angles input type
            Frac(f).PointM = Frac_Input(F).CenterCoord;  Frac(f).PointM(3)=0;

            % Obtaining the corners before rotation
            PointA_noRot       = [  Frac(f).PointM(1) - Frac_Input(F).Length /2  ;  Frac(f).PointM(2)  ;  0  ];
            PointB_noRot       = [  Frac(f).PointM(1) + Frac_Input(F).Length /2  ;  Frac(f).PointM(2)  ;  0  ];
            
            % Rotating around Z axis with horizontal rotation angle
            RotAxis_Hor        = [ 0 ; 0 ; 1 ];
            RotAngle_Hor_Rad   = Frac_Input(F).RotAngleAlongZ * pi / 180;
            Frac(f).PointA     = Rotate_Point_Around_Line_3D( PointA_noRot , RotAxis_Hor , Frac(f).PointM , RotAngle_Hor_Rad );
            Frac(f).PointB     = Rotate_Point_Around_Line_3D( PointB_noRot , RotAxis_Hor , Frac(f).PointM , RotAngle_Hor_Rad );

            % The coordinates of all Points in one variable (matrix) for Plotting
            Frac(f).Points = [ Frac(f).PointA , Frac(f).PointB , Frac(f).PointM ];

        end

        if ( min(Frac(f).Points(1,:)) < 0 ) || ( max(Frac(f).Points(1,:)) > LX ) || ...
           ( min(Frac(f).Points(2,:)) < 0 ) || ( max(Frac(f).Points(2,:)) > LY )
               error('Fracture #%1.0f does not fit in the matrix dimensions!',f);
        end

        %% Descretizing the Fracture Plates

        % The area of the fracture line (has two sides)
        Frac(f).Area_Total = norm( Frac(f).PointB - Frac(f).PointA ) * 2;

        % AB vector of the frwacture line
        Frac(f).AB_vec = Frac(f).PointB - Frac(f).PointA;

        % The equation of fracture line ax+by=c , using the central Point
        Frac(f).Equation.a =   Frac(f).AB_vec(2);
        Frac(f).Equation.b = - Frac(f).AB_vec(1);
        Frac(f).Equation.c =   Frac(f).AB_vec(2)*Frac(f).PointM(1) - Frac(f).AB_vec(1)*Frac(f).PointM(2);

        % Check if the fracture line is paralel to any of X,Y axes
        Frac(f).isParallel = string('None');
        if ( Frac(f).Equation.b == 0 ) , Frac(f).isParallel = string('AlongX');  end
        if ( Frac(f).Equation.a == 0 ) , Frac(f).isParallel = string('AlongY');  end

        % Length (Length_AB) of the fracture line
        Frac(f).Length_AB = norm( Frac(f).AB_vec );

        % Number of fractture line grid cells
        if ~isnan( Frac_Input(F).GridResRatio )
            Frac(f).N_Length_AB = round( ( Frac(f).Length_AB / max([LX,LY]) ) * max([NX,NY]) ^ Frac_Input(F).GridResRatio );
            Frac(f).N_Length_AB = max( Frac(f).N_Length_AB , 1);
        else
            Frac(f).N_Length_AB = Frac_Input(F).GridNumAlongL;
        end

        % The size of each fracture line grid cell with its vector
        Frac(f).D_Length_AB     = Frac(f).Length_AB / Frac(f).N_Length_AB;
        Frac(f).D_Length_AB_vec = Frac(f).AB_vec    / Frac(f).N_Length_AB;

        % Aperture, Porosity and Permeability of The fracture plate
        Frac(f).Aperture     = Frac_Input(F).Aperture;
        Frac(f).Porosity     = Frac_Input(F).Porosity;
        Frac(f).Permeability = Frac_Input(F).Permeability;

        % Coordinates (x,y,z) of cell centers and cell interfaces of the fracture line
        Frac(f).CellCenterCoords    = zeros( Frac(f).N_Length_AB   , 3 );
        Frac(f).GridCoords = zeros( Frac(f).N_Length_AB+1 , 3 );

        % Coordinates of fracture cell centers
        Frac(f).CellCenterCoords(:,1) = linspace( Frac(f).PointA(1) + (1/2)*Frac(f).D_Length_AB_vec(1) , Frac(f).PointB(1) - (1/2)*Frac(f).D_Length_AB_vec(1) , Frac(f).N_Length_AB );
        Frac(f).CellCenterCoords(:,2) = linspace( Frac(f).PointA(2) + (1/2)*Frac(f).D_Length_AB_vec(2) , Frac(f).PointB(2) - (1/2)*Frac(f).D_Length_AB_vec(2) , Frac(f).N_Length_AB );
        Frac(f).CellCenterCoords(:,3) = linspace( 0                                                    , 0                                                    , Frac(f).N_Length_AB );

        % Coordinates of fracture cell interfaces
        Frac(f).GridCoords(:,1) = linspace( Frac(f).PointA(1) , Frac(f).PointB(1) , Frac(f).N_Length_AB+1 );
        Frac(f).GridCoords(:,2) = linspace( Frac(f).PointA(2) , Frac(f).PointB(2) , Frac(f).N_Length_AB+1 );
        Frac(f).GridCoords(:,3) = linspace( 0                 , 0                 , Frac(f).N_Length_AB+1 );
        
    end
end

%% Intersections of fracture plate with each matrix cell (if any)
dummy = LX*10;                                                                        % A dummy value (an abnormal value)
for f = 1 : length(Frac)
        
    Frac(f).  intersectCoord_matCell = cell( Frac(f).N_Length_AB , 1 );    % Coordinates of intersections between each fracture cell and each matrix cell
    Frac(f).        areaFrac_matCell = cell( Frac(f).N_Length_AB , 1 );    % Area fraction of each fracture cell inside each matrix cube
    Frac(f).         aveDist_matCell = cell( Frac(f).N_Length_AB , 1 );    % Average distance between each fracture cell and each matrix cube
    Frac(f).    areaFrac_matCell_sum = 0;                                  % Summation of all area fractions (only to validate)
    Frac(f).  intersectCoord_fracObj = cell( length(Frac) , 1 );           % Coordinates of intersections between each two whole fracture plates
 
    Frac(f).    StarDelta2D_Neighbor = cell( length(Frac) , 1 );           % Indexes of neighboring fractures and thier cell in case of intersections
    
    for i_f = 1 : Frac(f).N_Length_AB
        Index_matIntersect = 0;

        % The matrix cell closest to fracture cell for the first intersection check
        [ ~ , i_1st ] = min( ( Frac(f).CellCenterCoords(i_f,1) - Xcm ).^2 );
        [ ~ , j_1st ] = min( ( Frac(f).CellCenterCoords(i_f,2) - Ycm ).^2 );
        
        Im_List  = Index_Matrix_2D( NX,NY , i_1st , j_1st );
        Im_count = 1;

        while Im_count <= length(Im_List)
            Im = Im_List(Im_count);
            Im_count = Im_count +1;
            
            % Retrieving i,j,k from the "Im" index             
            i = mod(  Im       , NX )   ;   if ( i==0 ),   i = NX;   end
            j = mod( (Im-i)/NX , NY ) +1;   if ( j==0 ),   j = NY;   end

            if Index_Matrix_2D(NX, NY, i, j) ~= Im
                error('i,j are not correspondent with Im. Check the formula again!');
            end
            
            doNotContinue        = 0;
            areCollinear         = 0;
            intersectCoord_temp  = [];
            intersectCoord_final = [];

            % Corner Points of the fracture line
            Line_fracCell.PointA = [ Frac(f).GridCoords(i_f   , 1) ; Frac(f).GridCoords(i_f   , 2) ; 0 ];
            Line_fracCell.PointB = [ Frac(f).GridCoords(i_f+1 , 1) ; Frac(f).GridCoords(i_f+1 , 2) ; 0 ];

            % Corner Points of the matrix Cell
            Plane_matCell.PointA = [ Xim(i  ) ; Yim(j  ) ; 0 ];
            Plane_matCell.PointB = [ Xim(i+1) ; Yim(j  ) ; 0 ];
            Plane_matCell.PointC = [ Xim(i+1) ; Yim(j+1) ; 0 ];
            Plane_matCell.PointD = [ Xim(i  ) ; Yim(j+1) ; 0 ];
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Matrix Cell, Line AB
            if doNotContinue == 0
                [lineGeostatus_AB, lineIntersectPoint_AB] = Line_Seg_Intersect_3D( Plane_matCell.PointA , Plane_matCell.PointB , Line_fracCell.PointA , Line_fracCell.PointB , almostZero );
                if lineGeostatus_AB.areCollinear == 1
                    doNotContinue = 1;
                    areCollinear  = 1;
                end
                intersectCoord_temp = [ intersectCoord_temp , lineIntersectPoint_AB ];
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Matrix Cell, Line BC
            if doNotContinue == 0
                [lineGeostatus_BC, lineIntersectPoint_BC] = Line_Seg_Intersect_3D( Plane_matCell.PointB , Plane_matCell.PointC , Line_fracCell.PointA , Line_fracCell.PointB , almostZero );
                if lineGeostatus_BC.areCollinear == 1
                    doNotContinue = 1;
                    areCollinear  = 1;
                end
                intersectCoord_temp = [ intersectCoord_temp , lineIntersectPoint_BC ];
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Matrix Cell, Line CD
            if doNotContinue == 0
                [lineGeostatus_CD, lineIntersectPoint_CD] = Line_Seg_Intersect_3D( Plane_matCell.PointC , Plane_matCell.PointD , Line_fracCell.PointA , Line_fracCell.PointB , almostZero );
                if lineGeostatus_CD.areCollinear == 1
                    doNotContinue = 1;
                    areCollinear  = 1;
                end
                intersectCoord_temp = [ intersectCoord_temp , lineIntersectPoint_CD ];
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Matrix Cell, Line DA
            if doNotContinue == 0
                [lineGeostatus_DA, lineIntersectPoint_DA] = Line_Seg_Intersect_3D( Plane_matCell.PointD , Plane_matCell.PointA , Line_fracCell.PointA , Line_fracCell.PointB , almostZero );
                if lineGeostatus_DA.areCollinear == 1
                    doNotContinue = 1;
                    areCollinear  = 1;
                end
                intersectCoord_temp = [ intersectCoord_temp , lineIntersectPoint_DA ];
            end
            
            intersectCoord_temp = [ intersectCoord_temp , Line_fracCell.PointA , Line_fracCell.PointB ];
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Processing the intersection Points
            for nr = 1 : size( intersectCoord_temp , 2 )

                % Removing the Points that are not inside the matrix cell
                if ( ( Xim(i) - intersectCoord_temp(1,nr) ) > almostZero ) || ( ( intersectCoord_temp(1,nr) - Xim(i+1) ) > almostZero )
                    intersectCoord_temp(:,nr) = [ dummy ; dummy ; dummy ];
                    continue;
                end
                if ( ( Yim(j) - intersectCoord_temp(2,nr) ) > almostZero ) || ( ( intersectCoord_temp(2,nr) - Yim(j+1) ) > almostZero )
                    intersectCoord_temp(:,nr) = [ dummy ; dummy ; dummy ];
                    continue;
                end

                % Removing the Points that are not inside the fracture cell
                if abs( Frac(f).Equation.a * intersectCoord_temp(1,nr) + Frac(f).Equation.b * intersectCoord_temp(2,nr) - Frac(f).Equation.c ) > almostZero
                    intersectCoord_temp(:,nr) = [ dummy ; dummy ; dummy ];
                    continue;
                end

                % Removing the Points that are repeated
                subtracted = ( intersectCoord_temp - intersectCoord_temp(:,nr) );
                subtracted = sqrt( subtracted(1,:).^2 + subtracted(2,:).^2 );
                if ~isempty( find( subtracted([1:nr-1 , nr+1:end]) < almostZero , 1 ) )
                    intersectCoord_temp(:,nr) = [ dummy ; dummy ; dummy ];
                    continue;
                end

                % Writing the intersection Points into a new variable
                if intersectCoord_temp(1,nr) ~= dummy
                    intersectCoord_final = [ intersectCoord_final , intersectCoord_temp(:,nr) ];
                end

            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Add the neighboring matrix cells to the list for intersection check if it is the first try
            % but no intersection ocuurs, so another one must be checked
            if ( ( isempty(intersectCoord_final) ) || ( size(intersectCoord_final,2) < 2 ) ) && Im_count == 2
                Im_next = Index_Matrix_2D( NX,NY,max(i-1,1 ),j           );  if ~ismember(Im_next,Im_List),  Im_List = [Im_List;Im_next];  end
                Im_next = Index_Matrix_2D( NX,NY,min(i+1,NX),j           );  if ~ismember(Im_next,Im_List),  Im_List = [Im_List;Im_next];  end
                Im_next = Index_Matrix_2D( NX,NY,i          ,max(j-1,1 ) );  if ~ismember(Im_next,Im_List),  Im_List = [Im_List;Im_next];  end
                Im_next = Index_Matrix_2D( NX,NY,i          ,min(j+1,NY) );  if ~ismember(Im_next,Im_List),  Im_List = [Im_List;Im_next];  end
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Assigning the index of matrix cells to the array
            if ( isempty(intersectCoord_final) ),  continue;  end

            Index_matIntersect = Index_matIntersect + 1;
            Frac(f).intersectCoord_matCell{i_f}{Index_matIntersect,1} = Im; 
            Frac(f).      areaFrac_matCell{i_f}{Index_matIntersect,1} = Im;
            Frac(f).       aveDist_matCell{i_f}{Index_matIntersect,1} = Im;

            % Add the neighboring cells to the list for intersection check
            Im_next = Index_Matrix_2D( NX,NY,max(i-1,1 ),j           );  if ~ismember(Im_next,Im_List),  Im_List = [Im_List;Im_next];  end
            Im_next = Index_Matrix_2D( NX,NY,min(i+1,NX),j           );  if ~ismember(Im_next,Im_List),  Im_List = [Im_List;Im_next];  end
            Im_next = Index_Matrix_2D( NX,NY,i          ,max(j-1,1 ) );  if ~ismember(Im_next,Im_List),  Im_List = [Im_List;Im_next];  end
            Im_next = Index_Matrix_2D( NX,NY,i          ,min(j+1,NY) );  if ~ismember(Im_next,Im_List),  Im_List = [Im_List;Im_next];  end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Writing the final intersection Points into the fracture data structure
            Frac(f).intersectCoord_matCell{i_f}{Index_matIntersect,2} = intersectCoord_final;
            if ( size(intersectCoord_final,2) > 2 )
                error('More than two intersection Points exist between the fracture cell and the matrix cell!');
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Calculating the area fraction
            Frac(f).areaFrac_matCell{i_f}{Index_matIntersect,2} = 0;
            if ( size(intersectCoord_final,2) < 2 )
                Frac(f).aveDist_matCell{i_f}{Index_matIntersect,2} = 1;  % If areaFrac_matCell is zero, make aveDist_matcube 1 to prevent 0/0 division.
                continue;
            end
            
            Frac(f).areaFrac_matCell{i_f}{Index_matIntersect,2} = norm( intersectCoord_final(:,1) - intersectCoord_final(:,2) );
            % Doubling the area fraction if the cell is not collinear to any matrix cube edges
            if areCollinear == 0
                Frac(f).areaFrac_matCell{i_f}{Index_matIntersect,2} = Frac(f).areaFrac_matCell{i_f}{Index_matIntersect,2} * 2;
            end
            Frac(f).areaFrac_matCell_sum = Frac(f).areaFrac_matCell_sum + Frac(f).areaFrac_matCell{i_f}{Index_matIntersect,2};

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Obtaining the average distance between the matrix cube cell and the fracture plate cell
            Frac(f).aveDist_matCell{i_f}{Index_matIntersect,2} = 0;
            [ ~ , Frac(f).aveDist_matCell{i_f}{Index_matIntersect,2} , ~] = Line_Plane_Connectivity_3D( Plane_matCell , intersectCoord_final(:,1) , intersectCoord_final(:,2) , almostZero );
           
        end % End of while-loop

    end % End of fractre 1st for-loop
    
end % End of main fracture for-loop

%% Plotting
for f = 1 : length(Frac)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plotting the boundaries of each fracture plate
    plot(Frac(f).Points(1,1:2),Frac(f).Points(2,1:2),'r','LineWidth',2);
    xlim([0,LX]); ylim([0,LY]);
    xlabel('X [cm]'); ylabel('Y [cm]');
    grid ON;
    hold on;
    
    % Writing the letters of A,B,C,D to the corners
    if length(Frac) <= 10
        text( Frac(f).PointA(1) , Frac(f).PointA(2) , 'A'        , 'Color' , 'black' , 'FontSize' , 14);
        text( Frac(f).PointB(1) , Frac(f).PointB(2) , 'B'        , 'Color' , 'black' , 'FontSize' , 14);
    end
    
    text( (2*Frac(f).PointA(1)+Frac(f).PointB(1))/3 , (2*Frac(f).PointA(2)+Frac(f).PointB(2))/3 , num2str(f) , 'Color' , 'black' , 'FontSize' , 14);
end

%% Intersections of fracture lines with each other (if any)
for f = 1 : length(Frac)
    for g = f+1 : length(Frac)
        

        % Obtaining the intersection Points
        [Geostatus_FF, intersectPoints_FF] = Line_Seg_Intersect_3D( Frac(f).PointA , Frac(f).PointB , Frac(g).PointA , Frac(g).PointB , almostZero );
        
        if ~isempty(intersectPoints_FF) && (Geostatus_FF.areCollinear ~= 1 )
            Frac(f).StarDelta2D_Neighbor{g} = struct;
            Frac(g).StarDelta2D_Neighbor{f} = struct;
        
            % Writing the final intersection Points into the fracture data structure
            Frac(f).intersectCoord_fracObj{g} = intersectPoints_FF;
            Frac(g).intersectCoord_fracObj{f} = intersectPoints_FF;
            
            distanceToRef = sqrt( sum( (Frac(f).GridCoords' - intersectPoints_FF).^2 , 1 ) );
            [ ~ , interface_index_f ] = min( distanceToRef );
            distanceToRef = sqrt( sum( (Frac(g).GridCoords' - intersectPoints_FF).^2 , 1 ) );
            [ ~ , interface_index_g ] = min( distanceToRef );
            
            if interface_index_f ~= 1
                Frac(f).StarDelta2D_Neighbor{g}.itself = [ interface_index_f-1 , interface_index_f ];
                Frac(g).StarDelta2D_Neighbor{f}.other  = [ interface_index_f-1 , interface_index_f ];
            else
                Frac(f).StarDelta2D_Neighbor{g}.itself = interface_index_f;
                Frac(g).StarDelta2D_Neighbor{f}.other  = interface_index_f;
            end
            if interface_index_g ~= 1
                Frac(g).StarDelta2D_Neighbor{f}.itself = [ interface_index_g-1 , interface_index_g ];
                Frac(f).StarDelta2D_Neighbor{g}.other  = [ interface_index_g-1 , interface_index_g ];
            else
                Frac(g).StarDelta2D_Neighbor{f}.itself = interface_index_g;
                Frac(f).StarDelta2D_Neighbor{g}.other  = interface_index_g;
            end

        end
        
    end
end

end