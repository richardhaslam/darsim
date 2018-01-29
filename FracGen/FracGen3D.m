%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FracGen3D  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Mousa HosseiniMehr, MSc Petroleum Engineering, CEG Faculty, TU Delft
% Project: 3D EDFM Package for F-ADM, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build Date : 2016-11-29
% Modified on: 2017-03-17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Frac = FracGen3D(obj)

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
Nt = 0;
%% Reporting ADM status of FracGen

%% Fractures Construction with Geometry Input
f = 0;
for F = 1 : size(Frac_Input,2)
    if Frac_Input(F).isActive == 1
            f = f + 1;

            if isnan( Frac_Input(F).CenterCoord(1) )
                % Fractures construction with three Points input type

                % The three vectors connecting 3 corners of the fracture plate (forming a triangle)
                Vec12 = Frac_Input(F).CornerCoords(:,2) - Frac_Input(F).CornerCoords(:,1);
                Vec23 = Frac_Input(F).CornerCoords(:,3) - Frac_Input(F).CornerCoords(:,2);
                Vec31 = Frac_Input(F).CornerCoords(:,3) - Frac_Input(F).CornerCoords(:,1);

                % Finding which of the three Points is the corner with the 90 deg angle              
                dotProduct = [ dot( Vec31 , Vec12 ) ; dot( Vec12 , Vec23 ) ; dot( Vec23 , Vec31 ) ];

                % Continue Here!!!
                PointRightAngle   = find( dotProduct==0 , 1 );
                PointSalientAngle = find( dotProduct~=0     );

                if isempty(PointRightAngle),  error('Fracture # %02d is not a rectangle! Please check the input file.',f);  end

                % Setting attributes to the corners (A is the corner with 90 deg angle)
                Frac(f).PointA   = Frac_Input(F).CornerCoords( : , PointRightAngle      );
                Frac(f).PointB   = Frac_Input(F).CornerCoords( : , PointSalientAngle(1) );
                Frac(f).PointD   = Frac_Input(F).CornerCoords( : , PointSalientAngle(2) );

                % The centeral Point of the plate 
                Frac(f).PointM = ( Frac(f).PointB + Frac(f).PointD ) /2;

                % Obtaining the 4th corner (by using the centeral Point)                         
                Frac(f).PointC = Frac(f).PointA + 2* ( Frac(f).PointM - Frac(f).PointA );

                % The coordinates of all Points in one variable (matrix) for Plotting
                Frac(f).Points = [ Frac(f).PointA , Frac(f).PointB , Frac(f).PointC , Frac(f).PointD , Frac(f).PointM ];

            else
                % Fracture construction with central Point and angles input type
                Frac(f).PointM = Frac_Input(F).CenterCoord;

                % Obtaining the corners before rotation
                PointA_noRot       = [  Frac(f).PointM(1) - Frac_Input(F).Length /2  ;  Frac(f).PointM(2) + Frac_Input(F).Width /2  ;  Frac(f).PointM(3)  ];
                PointB_noRot       = [  Frac(f).PointM(1) + Frac_Input(F).Length /2  ;  Frac(f).PointM(2) + Frac_Input(F).Width /2  ;  Frac(f).PointM(3)  ];
                PointC_noRot       = [  Frac(f).PointM(1) + Frac_Input(F).Length /2  ;  Frac(f).PointM(2) - Frac_Input(F).Width /2  ;  Frac(f).PointM(3)  ];
                PointD_noRot       = [  Frac(f).PointM(1) - Frac_Input(F).Length /2  ;  Frac(f).PointM(2) - Frac_Input(F).Width /2  ;  Frac(f).PointM(3)  ];   

                % Rotating around Z axis with horizontal rotation angle
                RotAxis_Hor        = [ 0 ; 0 ; 1 ];
                RotAngle_Hor_Rad   = Frac_Input(F).RotAngleAlongZ * pi / 180;
                PointA_Rot_Hor     = Rotate_Point_Around_Line_3D( PointA_noRot , RotAxis_Hor , Frac(f).PointM , RotAngle_Hor_Rad );
                PointB_Rot_Hor     = Rotate_Point_Around_Line_3D( PointB_noRot , RotAxis_Hor , Frac(f).PointM , RotAngle_Hor_Rad );
                PointC_Rot_Hor     = Rotate_Point_Around_Line_3D( PointC_noRot , RotAxis_Hor , Frac(f).PointM , RotAngle_Hor_Rad );
                PointD_Rot_Hor     = Rotate_Point_Around_Line_3D( PointD_noRot , RotAxis_Hor , Frac(f).PointM , RotAngle_Hor_Rad );

                % Rotating around "M_AB - M_CD" with along width rotation angle
                PointBMC           = ( PointB_Rot_Hor + PointC_Rot_Hor ) /2;
                PointDMA           = ( PointD_Rot_Hor + PointA_Rot_Hor ) /2;
                RotAxis_Length     = PointDMA - PointBMC;
                RotAxis_Length_Rad  = Frac_Input(F).RotAngleAlongW * pi / 180;
                PointA_Rot_Length   = Rotate_Point_Around_Line_3D( PointA_Rot_Hor , RotAxis_Length , Frac(f).PointM , RotAxis_Length_Rad );
                PointB_Rot_Length   = Rotate_Point_Around_Line_3D( PointB_Rot_Hor , RotAxis_Length , Frac(f).PointM , RotAxis_Length_Rad );
                PointC_Rot_Length   = Rotate_Point_Around_Line_3D( PointC_Rot_Hor , RotAxis_Length , Frac(f).PointM , RotAxis_Length_Rad );
                PointD_Rot_Length   = Rotate_Point_Around_Line_3D( PointD_Rot_Hor , RotAxis_Length , Frac(f).PointM , RotAxis_Length_Rad );

                % Rotating around "M_BC - M_DA" with along length rotation angle
                PointAMB           = ( PointA_Rot_Length + PointB_Rot_Length ) /2;
                PointCMD           = ( PointC_Rot_Length + PointD_Rot_Length ) /2;
                RotAxis_Width      = PointCMD - PointAMB;
                RotAxis_Width_Rad = Frac_Input(F).RotAngleAlongL * pi / 180;
                Frac(f).PointA     = Rotate_Point_Around_Line_3D( PointA_Rot_Length , RotAxis_Width , Frac(f).PointM , RotAxis_Width_Rad );
                Frac(f).PointB     = Rotate_Point_Around_Line_3D( PointB_Rot_Length , RotAxis_Width , Frac(f).PointM , RotAxis_Width_Rad );
                Frac(f).PointC     = Rotate_Point_Around_Line_3D( PointC_Rot_Length , RotAxis_Width , Frac(f).PointM , RotAxis_Width_Rad );
                Frac(f).PointD     = Rotate_Point_Around_Line_3D( PointD_Rot_Length , RotAxis_Width , Frac(f).PointM , RotAxis_Width_Rad );

                % The coordinates of all Points in one variable (matrix) for Plotting
                Frac(f).Points = [ Frac(f).PointA , Frac(f).PointB , Frac(f).PointC , Frac(f).PointD , Frac(f).PointM ];

            end

            if ( min(Frac(f).Points(1,:)) < 0 - almostZero ) || ( max(Frac(f).Points(1,:)) > LX + almostZero ) || ...
               ( min(Frac(f).Points(2,:)) < 0 - almostZero ) || ( max(Frac(f).Points(2,:)) > LY + almostZero ) || ...
               ( min(Frac(f).Points(3,:)) < 0 - almostZero ) || ( max(Frac(f).Points(3,:)) > LZ + almostZero )
                   error('Fracture #%02d does not fit in the matrix dimensions! Please check the input file.',f);
            end

    %% Descretizing the Fracture Plates

        % The area of the fracture plate (each plate has two sides)
        Frac(f).Area_Total = Triangle_Area_3D( Frac(f).PointA , Frac(f).PointB , Frac(f).PointC ) * 2 * 2;

        % AB and AD vectors of the plate
        Frac(f).AB_vec = Frac(f).PointB - Frac(f).PointA;
        Frac(f).AD_vec = Frac(f).PointD - Frac(f).PointA;
        Frac(f).AC_vec = Frac(f).PointC - Frac(f).PointA;

        % The normal vector of plate
        Frac(f).n_vec = cross( Frac(f).AB_vec , Frac(f).AD_vec );
        Frac(f).n_vec = Frac(f).n_vec / norm( Frac(f).n_vec );

        % The normal vector of plate with its size normalized to 1/10 of the reservoir length (LX) (just for plotting)
        Frac(f).n_vec_plot = Frac(f).n_vec * ( Frac_Input(F).Length + Frac_Input(F).Width ) /10;

        % The equation of fracture plate a(x-x_M)+b(y-y_M)+c(z-z_M)=0 or ax+by+cz=d , using the central Point
        Frac(f).Equation.a = Frac(f).n_vec(1);
        Frac(f).Equation.b = Frac(f).n_vec(2);
        Frac(f).Equation.c = Frac(f).n_vec(3);
        Frac(f).Equation.d = Frac(f).Equation.a*Frac(f).PointM(1) + Frac(f).Equation.b*Frac(f).PointM(2) + Frac(f).Equation.c*Frac(f).PointM(3);

        % Check if the fracture plate is paralel to any of X,Y,Z axes
        Frac(f).isParallel = string('None');
        if ( Frac(f).Equation.b == 0 ) && ( Frac(f).Equation.c == 0), Frac(f).isParallel = string('AlongX');  end
        if ( Frac(f).Equation.a == 0 ) && ( Frac(f).Equation.c == 0), Frac(f).isParallel = string('AlongY');  end
        if ( Frac(f).Equation.a == 0 ) && ( Frac(f).Equation.b == 0), Frac(f).isParallel = string('AlongZ');  end

        % Length (Length_AB) and width (Width_AD) of the fracture plate
        Frac(f).Length_AB = norm( Frac(f).AB_vec );
        Frac(f).Width_AD  = norm( Frac(f).AD_vec );

        % Adding ADM configuration
        Frac(f).ADM = Frac_Input(F).ADM;
        
        % Number of fractture plate grid cells in AB and AD directions
        if ~isnan( Frac_Input(F).GridResRatio )
            Frac(f).N_Length_AB = round( ( Frac(f).Length_AB / max([LX,LY,LZ]) ) * max([NX,NY,NZ]) ^ Frac_Input(F).GridResRatio );
            Frac(f).N_Width_AD  = round( ( Frac(f).Width_AD  / max([LX,LY,LZ]) ) * max([NX,NY,NZ]) ^ Frac_Input(F).GridResRatio );
            if Frac(f).ADM(1)
                if ( mod( Frac(f).N_Length_AB , Frac(f).ADM(3)^Frac(f).ADM(2) )~=0 ) % Correcting N_Length_AB to match MMs level and coarsening ratio
                    Frac(f).N_Length_AB = Frac(f).N_Length_AB + Frac(f).ADM(3)^Frac(f).ADM(2) - mod( Frac(f).N_Length_AB , Frac(f).ADM(3)^Frac(f).ADM(2) );
                end
                Frac(f).N_Length_AB = max( Frac(f).N_Length_AB , 3*Frac(f).ADM(3)^Frac(f).ADM(2) ); % At least 3 Coarse Cells
            end
            if Frac(f).ADM(1)
                if ( mod( Frac(f).N_Width_AD , Frac(f).ADM(4)^Frac(f).ADM(2) )~=0 ) % Correcting N_Width_AD to match MMs level and coarsening ratio
                    Frac(f).N_Width_AD = Frac(f).N_Width_AD + Frac(f).ADM(4)^Frac(f).ADM(2) - mod( Frac(f).N_Width_AD , Frac(f).ADM(4)^Frac(f).ADM(2) );
                end
                Frac(f).N_Width_AD = max( Frac(f).N_Width_AD , 3*Frac(f).ADM(4)^Frac(f).ADM(2) ); % At least 3 Coarse Cells
            end
            if NZ == 1,  Frac(f).N_Width_AD = 1;  end
            Frac(f).N_Length_AB = max( Frac(f).N_Length_AB , 1);
            Frac(f).N_Width_AD  = max( Frac(f).N_Width_AD  , 1);
        else
            Frac(f).N_Length_AB = Frac_Input(F).GridNumAlongL;
            Frac(f).N_Width_AD  = Frac_Input(F).GridNumAlongW;
            if Frac(f).ADM(1)
                if mod( Frac(f).N_Length_AB , Frac(f).ADM(3)^Frac(f).ADM(2) ) ~= 0
                    error('In fracture #%02d, ADM coarsening ratio along length is not acceptable to the number of grid cells! Please check the input file.',f);
                end
                if mod( Frac(f).N_Width_AD  , Frac(f).ADM(4)^Frac(f).ADM(2) ) ~= 0
                    error('In fracture #%02d, ADM coarsening ratio along width is not acceptable to the number of grid cells! Please check the input file.',f);
                end
            end
        end
        
        % The size of each fracture plate grid cell in AB and AD directions with their vectors
        Frac(f).D_Length_AB     = Frac(f).Length_AB / Frac(f).N_Length_AB;
        Frac(f).D_Width_AD     = Frac(f).Width_AD  / Frac(f).N_Width_AD;
        Frac(f).D_Length_AB_vec = Frac(f).AB_vec    / Frac(f).N_Length_AB;
        Frac(f).D_Width_AD_vec = Frac(f).AD_vec    / Frac(f).N_Width_AD;

        % Aperture, Porosity and Permeability of The fracture plate
        Frac(f).Aperture     = Frac_Input(F).Aperture;
        Frac(f).Porosity     = Frac_Input(F).Porosity;
        Frac(f).Permeability = Frac_Input(F).Permeability;

        % Coordinates (x,y,z) of cell centers and cell interfaces of the fracture plate
        Frac(f).CellCenterCoords    = zeros( Frac(f).N_Length_AB   , Frac(f).N_Width_AD   , 3 );
        Frac(f).CellCenterCoordsV2  = zeros( Frac(f).N_Length_AB * Frac(f).N_Width_AD , 3 );
        Frac(f).GridCoords = zeros( Frac(f).N_Length_AB+1 , Frac(f).N_Width_AD+1 , 3 );

        % Coordinates of fracture cell centers
        for i = 1 : Frac(f).N_Length_AB
            Frac(f).CellCenterCoords(i,:,1) = linspace( Frac(f).PointA(1) + (i-1/2)*Frac(f).D_Length_AB_vec(1) + Frac(f).D_Width_AD_vec(1)/2 , ...
                                                   Frac(f).PointD(1) + (i-1/2)*Frac(f).D_Length_AB_vec(1) - Frac(f).D_Width_AD_vec(1)/2 , Frac(f).N_Width_AD );
            Frac(f).CellCenterCoords(i,:,2) = linspace( Frac(f).PointA(2) + (i-1/2)*Frac(f).D_Length_AB_vec(2) + Frac(f).D_Width_AD_vec(2)/2 , ...
                                                   Frac(f).PointD(2) + (i-1/2)*Frac(f).D_Length_AB_vec(2) - Frac(f).D_Width_AD_vec(2)/2 , Frac(f).N_Width_AD );
            Frac(f).CellCenterCoords(i,:,3) = linspace( Frac(f).PointA(3) + (i-1/2)*Frac(f).D_Length_AB_vec(3) + Frac(f).D_Width_AD_vec(3)/2 , ...
                                                   Frac(f).PointD(3) + (i-1/2)*Frac(f).D_Length_AB_vec(3) - Frac(f).D_Width_AD_vec(3)/2 , Frac(f).N_Width_AD );
        end
        
        % Coordinates of fracture cell centers in vector form
        for j = 1 : Frac(f).N_Width_AD
            Frac(f).CellCenterCoordsV2( (j-1)*Frac(f).N_Length_AB+1 : j*Frac(f).N_Length_AB , : ) = Frac(f).CellCenterCoords(:,j,:);
        end

        % Coordinates of fracture cell interfaces
        for i = 1 : Frac(f).N_Length_AB+1
            Frac(f).GridCoords(i,:,1) = linspace( Frac(f).PointA(1) + (i-1)*Frac(f).D_Length_AB_vec(1) , ...
                                                      Frac(f).PointD(1) + (i-1)*Frac(f).D_Length_AB_vec(1) , Frac(f).N_Width_AD+1 );
            Frac(f).GridCoords(i,:,2) = linspace( Frac(f).PointA(2) + (i-1)*Frac(f).D_Length_AB_vec(2) , ...
                                                      Frac(f).PointD(2) + (i-1)*Frac(f).D_Length_AB_vec(2) , Frac(f).N_Width_AD+1 );
            Frac(f).GridCoords(i,:,3) = linspace( Frac(f).PointA(3) + (i-1)*Frac(f).D_Length_AB_vec(3) , ...
                                                      Frac(f).PointD(3) + (i-1)*Frac(f).D_Length_AB_vec(3) , Frac(f).N_Width_AD+1 );    
        end
        
        % Reporting fracture properties:
        fprintf('Fracture %2d: Dimension= %5.2f x %5.2f [m2] , Grid= %3.0f x %3.0f = %4.0f , ADM lvl= %1.0f\n', ...
        f, Frac(f).Length_AB, Frac(f).Width_AD, Frac(f).N_Length_AB, Frac(f).N_Width_AD, Frac(f).N_Length_AB*Frac(f).N_Width_AD, Frac(f).ADM(1)*Frac(f).ADM(2) );
        
        Nt = Nt + Frac(f).N_Length_AB*Frac(f).N_Width_AD;
    end
end
% Row-wise to column-wise
Frac = Frac';
%
fprintf('\n');
fprintf('Total number of fractures: %3.0f\n',length(Frac));
fprintf('Total number of fractures cells: %3.0f\n',Nt);
fprintf('---------------------------------------------------------\n');

%% Plotting fracture plates
fprintf('Plotting fractures ...\n');
figure(1); 
for f = 1 : length(Frac)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plotting the boundaries of each fracture plate
    plot3(Frac(f).Points(1,[1:4,1]),Frac(f).Points(2,[1:4,1]),Frac(f).Points(3,[1:4,1]),'r','LineWidth',4);
    xlim([0,LX]); ylim([0,LY]); zlim([0,LZ]);
    xlabel('X [cm]'); ylabel('Y [cm]'); zlabel('Z [cm]');
    grid ON;
    
    % Plotting the normal vector of each fracture plate
    hold on;
    quiver3(Frac(f).PointM(1),Frac(f).PointM(2),Frac(f).PointM(3),Frac(f).n_vec_plot(1),Frac(f).n_vec_plot(2),Frac(f).n_vec_plot(3),'b','LineWidth',3);
    
    % Plotting cell edges of each fracture plate
    for i = 1 : Frac(f).N_Length_AB+1
        plot3(Frac(f).GridCoords(i,:,1),Frac(f).GridCoords(i,:,2),Frac(f).GridCoords(i,:,3),'r','LineWidth',2);
    end
    for j = 1 : Frac(f).N_Width_AD+1
        plot3(Frac(f).GridCoords(:,j,1),Frac(f).GridCoords(:,j,2),Frac(f).GridCoords(:,j,3),'r','LineWidth',2);
    end
    
    % Writing the letters of A,B,C,D to the corners
    text( Frac(f).PointA(1) - 0.05*Frac(f).AD_vec(1) - 0.05*Frac(f).AB_vec(1) , ...
          Frac(f).PointA(2) - 0.05*Frac(f).AD_vec(2) - 0.05*Frac(f).AB_vec(2) , ...
          Frac(f).PointA(3) - 0.05*Frac(f).AD_vec(3) - 0.05*Frac(f).AB_vec(3) , ...
          'A' , 'Color' , 'black' , 'FontSize' , 14);
    text( Frac(f).PointB(1) - 0.05*Frac(f).AD_vec(1) + 0.05*Frac(f).AB_vec(1) , ...
          Frac(f).PointB(2) - 0.05*Frac(f).AD_vec(2) + 0.05*Frac(f).AB_vec(2) , ...
          Frac(f).PointB(3) - 0.05*Frac(f).AD_vec(3) + 0.05*Frac(f).AB_vec(3) , ...
          'B' , 'Color' , 'black' , 'FontSize' , 14);
    text( Frac(f).PointC(1) + 0.05*Frac(f).AD_vec(1) + 0.05*Frac(f).AB_vec(1) , ...
          Frac(f).PointC(2) + 0.05*Frac(f).AD_vec(2) + 0.05*Frac(f).AB_vec(2) , ...
          Frac(f).PointC(3) + 0.05*Frac(f).AD_vec(3) + 0.05*Frac(f).AB_vec(3) , ...
          'C' , 'Color' , 'black' , 'FontSize' , 14);
    text( Frac(f).PointD(1) + 0.05*Frac(f).AD_vec(1) - 0.05*Frac(f).AB_vec(1) , ...
          Frac(f).PointD(2) + 0.05*Frac(f).AD_vec(2) - 0.05*Frac(f).AB_vec(2) , ...
          Frac(f).PointD(3) + 0.05*Frac(f).AD_vec(3) - 0.05*Frac(f).AB_vec(3) , ...
          'D' , 'Color' , 'black' , 'FontSize' , 14);
end

% Plotting the edges of matrix cells
if max([NX,NY,NZ]) <= 10
    hold on;
    edgeColor = [ 0 0 0 0.25];
    for i = 1:NX+1
        for j = 1:NY+1
            plot3( [Xim(i),Xim(i)] , [Yim(j),Yim(j)] , [Zim(1),Zim(end)] , 'Color' , edgeColor );
        end
    end
    for i = 1:NX+1
        for k = 1:NZ+1
            plot3( [Xim(i),Xim(i)] , [Yim(1),Yim(end)] , [Zim(k),Zim(k)] , 'Color' , edgeColor );
        end
    end
    for j = 1:NY+1
        for k = 1:NZ+1
            plot3( [Xim(1),Xim(end)] , [Yim(j),Yim(j)] , [Zim(k),Zim(k)] , 'Color' , edgeColor );
        end
    end
end

rotate3d on;
view([0 90]);
fprintf('---------------------------------------------------------\n');

%% Intializing
for f = 1 : length(Frac)
    Frac(f).  intersectCoord_fracObj = cell( length(Frac) , 1 );           % Coordinates of intersections between each two whole fracture plates
    Frac(f).            Overlap_frac = cell( length(Frac) , 1 );           % Overlap status of fracture f cells due to intersection with fracture g
end
%return;
%% Intersections of fracture plates with each other (if any)
fprintf('Obtaining the intersection of fracture plates ...\n');
dummy = LX*10;                                                             % A dummy value (an abnormal value)
for f = 1 : length(Frac)
    for g = f+1 : length(Frac)
        
        % Obtaining the intersection Points
        [Geostatus_PP, intersectPoints] = Plane_Seg_Intersect_3D( Frac(f) , Frac(g) , almostZero );
        intersectCoord_fracObj_final = [];
        
        if ~isempty(intersectPoints) && (Geostatus_PP.areCoplanar ~= 1 )
            % Processing the intersection Points
            for nr = 1 : size(intersectPoints,2)
                
                % Removing the Points that are not inside the fracture rectangles
                isInside1 = Is_Point_Inside_Rectangle_3D( Frac(f) , intersectPoints(:,nr) , almostZero );
                isInside2 = Is_Point_Inside_Rectangle_3D( Frac(g) , intersectPoints(:,nr) , almostZero );
                if ( isInside1 == 0 ) || ( isInside2 == 0 )
                    intersectPoints(:,nr) = [ dummy ; dummy ; dummy ];
                    continue;
                end
                
                % Removing the points that are repeated
                subtracted = ( intersectPoints - intersectPoints(:,nr) );
                subtracted = sqrt( subtracted(1,:).^2 + subtracted(2,:).^2 + subtracted(3,:).^2 );
                if ~isempty( find( subtracted([1:nr-1 , nr+1:end]) <= almostZero , 1 ) )
                    intersectPoints(:,nr) = [ dummy ; dummy ; dummy ];
                    continue;
                end
                
                if intersectPoints(1,nr) ~= dummy
                    intersectCoord_fracObj_final = [ intersectCoord_fracObj_final , intersectPoints(:,nr) ];
                end
            end
            
            % If there is no or only one final intersection Point, do not plot
            if size(intersectCoord_fracObj_final,2) < 2,  continue;  end
            
            % Writing the final intersection Points into the fracture data structure
            Frac(f).intersectCoord_fracObj{g} = intersectCoord_fracObj_final;
            Frac(g).intersectCoord_fracObj{f} = intersectCoord_fracObj_final;
            Frac(f).          Overlap_frac{g} = zeros( Frac(f).N_Length_AB * Frac(f).N_Width_AD , 1 );
            Frac(g).          Overlap_frac{f} = zeros( Frac(g).N_Length_AB * Frac(g).N_Width_AD , 1 );
            
            % Plotting the intersection line of each two fracture plates
            figure(1);
            plot3(intersectCoord_fracObj_final(1,:),intersectCoord_fracObj_final(2,:),intersectCoord_fracObj_final(3,:),'Color',[0,0.5,0],'LineWidth',5);
            
        end
    end
end
fprintf('---------------------------------------------------------\n');

%% Intersections of fracture plate with each matrix cell
fprintf('Obtaining fractures - matrix overlaps:\n',f);
fprintf('---> Fracture ');
for f = 1 : length(Frac)
    if (f>1),  fprintf(repmat('\b', 1, 6));  end
    fprintf('%02d/%02d\n',f,length(Frac));
    Frac(f).  intersectCoord_matCell = cell( Frac(f).N_Length_AB*Frac(f).N_Width_AD , 1 );   % Coordinates of intersections between each fracture cell and each matrix cell
    Frac(f).        areaFrac_matCell = cell( Frac(f).N_Length_AB*Frac(f).N_Width_AD , 1 );   % Area fraction of each fracture cell inside each matrix cube
    Frac(f).         aveDist_matCell = cell( Frac(f).N_Length_AB*Frac(f).N_Width_AD , 1 );   % Average distance between each fracture cell and each matrix cube
    Frac(f).           T_Geo_matCell = cell( Frac(f).N_Length_AB*Frac(f).N_Width_AD , 1 );   % Geometrical Transmissibility between each fracture cell and each matrix cube
    Frac(f).    areaFrac_matCell_sum = 0;                                                    % Summation of all area fractions (only to validate)
    Frac(f). intersectCoord_fracCell = cell( length(Frac) , 1 );                             % Coordinates of intersections between each two fracture cells of distinct fracture plates
    Frac(f).       areaFrac_fracCell = cell( length(Frac) , 1 );                             % Area fraction of of intersection line between each two fracture cells of distinct fracture plates
    Frac(f).        aveDist_fracCell = cell( length(Frac) , 1 );                             % Average distance between each fracture cell and the intersection line between each two fracture cells of distinct fracture plates
    Frac(f).          T_Geo_fracCell = cell( length(Frac) , 1 );                             % Geometrical Transmissibility between each two non-neighboring fracture cells of distinct fracture plates
    Frac(f).      NumOf_fracCellConn = zeros( Frac(f).N_Length_AB*Frac(f).N_Width_AD , 1);   % Number of non-neighboring fracture cells connections
    
    for j_f = 1 : Frac(f).N_Width_AD
        for i_f = 1 : Frac(f).N_Length_AB
            
            If = Index_Fracture_2D( Frac(f).N_Length_AB , Frac(f).N_Width_AD , i_f , j_f );
            Index_matIntersect = 0;
            
            % The matrix cell closest to fracture cell for the first intersection check
            [ ~ , i_1st ] = min( ( Frac(f).CellCenterCoords(i_f,j_f,1) - Xcm ).^2 );
            [ ~ , j_1st ] = min( ( Frac(f).CellCenterCoords(i_f,j_f,2) - Ycm ).^2 );
            [ ~ , k_1st ] = min( ( Frac(f).CellCenterCoords(i_f,j_f,3) - Zcm ).^2 );
            Im_List  = Index_Matrix_3D( NX,NY,NZ , i_1st , j_1st , k_1st );
            Im_count = 1;
            
            while Im_count <= length(Im_List)
                Im = Im_List(Im_count);
                Im_count = Im_count +1;
                
                % Retrieving i,j,k from the "Im" index             
                i = mod(   Im                , NX )   ;   if ( i==0 ),   i = NX;   end
                j = mod(  (Im-i)/NX          , NY ) +1;   if ( j==0 ),   j = NY;   end
                k = mod( ((Im-i)/NX -j+1)/NY , NZ ) +1;   if ( k==0 ),   k = NZ;   end
                if Index_Matrix_3D(NX, NY, NZ, i, j, k) ~= Im
                    error('i,j,k are not correspondent with Im. Check the formula again!');
                end

                doNotContinue        = 0;
                areCoplanar          = 0;
                intersectCoord_temp  = [];
                intersectCoord_final = [];

                % Corner Points of the fracture cell
                Plane_fracCell.PointA = [ Frac(f).GridCoords(i_f  ,j_f  ,1) ; Frac(f).GridCoords(i_f  ,j_f  ,2) ; Frac(f).GridCoords(i_f  ,j_f  ,3) ];
                Plane_fracCell.PointB = [ Frac(f).GridCoords(i_f  ,j_f+1,1) ; Frac(f).GridCoords(i_f  ,j_f+1,2) ; Frac(f).GridCoords(i_f  ,j_f+1,3) ];
                Plane_fracCell.PointC = [ Frac(f).GridCoords(i_f+1,j_f+1,1) ; Frac(f).GridCoords(i_f+1,j_f+1,2) ; Frac(f).GridCoords(i_f+1,j_f+1,3) ];
                Plane_fracCell.PointD = [ Frac(f).GridCoords(i_f+1,j_f  ,1) ; Frac(f).GridCoords(i_f+1,j_f  ,2) ; Frac(f).GridCoords(i_f+1,j_f  ,3) ];

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Matrix Cube Face #X1
                if doNotContinue == 0
                    Plane_matFace_X1.PointA = [ Xim(i  ) ; Yim(j  ) ; Zim(k  ) ];
                    Plane_matFace_X1.PointB = [ Xim(i  ) ; Yim(j+1) ; Zim(k  ) ];
                    Plane_matFace_X1.PointC = [ Xim(i  ) ; Yim(j+1) ; Zim(k+1) ];
                    Plane_matFace_X1.PointD = [ Xim(i  ) ; Yim(j  ) ; Zim(k+1) ];
                    
                    [Geostatus_X1, intersectPoints_X1] = Plane_Seg_Intersect_3D( Plane_matFace_X1 , Plane_fracCell , almostZero );
                    if ( Geostatus_X1.areCoplanar == 1 )
                        doNotContinue      = 1;
                        areCoplanar        = 1;
                        intersectPoints_X1 = [ intersectPoints_X1 , Plane_matFace_X1.PointA , Plane_matFace_X1.PointB , Plane_matFace_X1.PointC , Plane_matFace_X1.PointD ];
                    end
                    intersectCoord_temp = [ intersectCoord_temp , intersectPoints_X1 ];
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Matrix Cube Face #X2
                if doNotContinue == 0
                    Plane_matFace_X2.PointA = [ Xim(i+1) ; Yim(j  ) ; Zim(k  ) ];
                    Plane_matFace_X2.PointB = [ Xim(i+1) ; Yim(j+1) ; Zim(k  ) ];
                    Plane_matFace_X2.PointC = [ Xim(i+1) ; Yim(j+1) ; Zim(k+1) ];
                    Plane_matFace_X2.PointD = [ Xim(i+1) ; Yim(j  ) ; Zim(k+1) ];

                    [Geostatus_X2, intersectPoints_X2] = Plane_Seg_Intersect_3D( Plane_matFace_X2 , Plane_fracCell , almostZero );
                    if ( Geostatus_X2.areCoplanar == 1 )
                        doNotContinue      = 1;
                        areCoplanar        = 1;
                        intersectPoints_X2 = [ intersectPoints_X2 , Plane_matFace_X2.PointA , Plane_matFace_X2.PointB , Plane_matFace_X2.PointC , Plane_matFace_X2.PointD ];
                    end
                    intersectCoord_temp = [ intersectCoord_temp , intersectPoints_X2 ];
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Matrix Cube Face #Y1
                if doNotContinue == 0
                    Plane_matFace_Y1.PointA = [ Xim(i  ) ; Yim(j  ) ; Zim(k  ) ];
                    Plane_matFace_Y1.PointB = [ Xim(i+1) ; Yim(j  ) ; Zim(k  ) ];
                    Plane_matFace_Y1.PointC = [ Xim(i+1) ; Yim(j  ) ; Zim(k+1) ];
                    Plane_matFace_Y1.PointD = [ Xim(i  ) ; Yim(j  ) ; Zim(k+1) ];

                    [Geostatus_Y1, intersectPoints_Y1] = Plane_Seg_Intersect_3D( Plane_matFace_Y1 , Plane_fracCell , almostZero );
                    if ( Geostatus_Y1.areCoplanar == 1 )
                        doNotContinue      = 1;
                        areCoplanar        = 1;
                        intersectPoints_Y1 = [ intersectPoints_Y1 , Plane_matFace_Y1.PointA , Plane_matFace_Y1.PointB , Plane_matFace_Y1.PointC , Plane_matFace_Y1.PointD ];
                    end
                    intersectCoord_temp = [ intersectCoord_temp , intersectPoints_Y1 ];
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Matrix Cube Face #Y2
                if doNotContinue == 0
                    Plane_matFace_Y2.PointA = [ Xim(i  ) ; Yim(j+1) ; Zim(k  ) ];
                    Plane_matFace_Y2.PointB = [ Xim(i+1) ; Yim(j+1) ; Zim(k  ) ];
                    Plane_matFace_Y2.PointC = [ Xim(i+1) ; Yim(j+1) ; Zim(k+1) ];
                    Plane_matFace_Y2.PointD = [ Xim(i  ) ; Yim(j+1) ; Zim(k+1) ];

                    [Geostatus_Y2, intersectPoints_Y2] = Plane_Seg_Intersect_3D( Plane_matFace_Y2 , Plane_fracCell , almostZero );
                    if ( Geostatus_Y2.areCoplanar == 1 )
                        doNotContinue      = 1;
                        areCoplanar        = 1;
                        intersectPoints_Y2 = [ intersectPoints_Y2 , Plane_matFace_Y2.PointA , Plane_matFace_Y2.PointB , Plane_matFace_Y2.PointC , Plane_matFace_Y2.PointD ];
                    end
                    intersectCoord_temp = [ intersectCoord_temp , intersectPoints_Y2 ];
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Matrix Cube Face #Z1
                if doNotContinue == 0
                    Plane_matFace_Z1.PointA = [ Xim(i  ) ; Yim(j  ) ; Zim(k  ) ];
                    Plane_matFace_Z1.PointB = [ Xim(i+1) ; Yim(j  ) ; Zim(k  ) ];
                    Plane_matFace_Z1.PointC = [ Xim(i+1) ; Yim(j+1) ; Zim(k  ) ];
                    Plane_matFace_Z1.PointD = [ Xim(i  ) ; Yim(j+1) ; Zim(k  ) ];

                    [Geostatus_Z1, intersectPoints_Z1] = Plane_Seg_Intersect_3D( Plane_matFace_Z1 , Plane_fracCell , almostZero );
                    if ( Geostatus_Z1.areCoplanar == 1 )
                        doNotContinue      = 1;
                        areCoplanar        = 1;
                        intersectPoints_Z1 = [ intersectPoints_Z1 , Plane_matFace_Z1.PointA , Plane_matFace_Z1.PointB , Plane_matFace_Z1.PointC , Plane_matFace_Z1.PointD ];
                    end
                    intersectCoord_temp = [ intersectCoord_temp , intersectPoints_Z1 ];
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Matrix Cube Face #Z2
                if doNotContinue == 0
                    Plane_matFace_Z2.PointA = [ Xim(i  ) ; Yim(j  ) ; Zim(k+1) ];
                    Plane_matFace_Z2.PointB = [ Xim(i+1) ; Yim(j  ) ; Zim(k+1) ];
                    Plane_matFace_Z2.PointC = [ Xim(i+1) ; Yim(j+1) ; Zim(k+1) ];
                    Plane_matFace_Z2.PointD = [ Xim(i  ) ; Yim(j+1) ; Zim(k+1) ];

                    [Geostatus_Z2, intersectPoints_Z2] = Plane_Seg_Intersect_3D( Plane_matFace_Z2 , Plane_fracCell , almostZero );
                    if ( Geostatus_Z2.areCoplanar == 1 )
                        doNotContinue = 1;
                        areCoplanar   = 1;
                        intersectPoints_Z2 = [ intersectPoints_Z2 , Plane_matFace_Z2.PointA , Plane_matFace_Z2.PointB , Plane_matFace_Z2.PointC , Plane_matFace_Z2.PointD ];
                    end
                    intersectCoord_temp = [ intersectCoord_temp , intersectPoints_Z2 ];
                end

                intersectCoord_temp = [ intersectCoord_temp , Plane_fracCell.PointA , Plane_fracCell.PointB , Plane_fracCell.PointC , Plane_fracCell.PointD ];

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
                    if ( ( Zim(k) - intersectCoord_temp(3,nr) ) > almostZero ) || ( ( intersectCoord_temp(3,nr) - Zim(k+1) ) > almostZero )
                        intersectCoord_temp(:,nr) = [ dummy ; dummy ; dummy ];
                        continue;
                    end

                    % Removing the Points that are not inside the fracture cell
                    isInside = Is_Point_Inside_Rectangle_3D( Plane_fracCell , intersectCoord_temp(:,nr) , almostZero );
                    if isInside == 0
                        intersectCoord_temp(:,nr) = [ dummy ; dummy ; dummy ];
                        continue;
                    end

                    % Removing the Points that are repeated
                    subtracted = ( intersectCoord_temp - intersectCoord_temp(:,nr) );
                    subtracted = sqrt( subtracted(1,:).^2 + subtracted(2,:).^2 + subtracted(3,:).^2 );
                    if ~isempty( find( subtracted([1:nr-1 , nr+1:end]) <= almostZero , 1 ) )
                        intersectCoord_temp(:,nr) = [ dummy ; dummy ; dummy ];
                        continue;
                    end

                    % Writing the intersection Points into a new variable
                    if intersectCoord_temp(1,nr) ~= dummy
                        intersectCoord_final = [ intersectCoord_final , intersectCoord_temp(:,nr) ];
                    end

                end
                intersectCoord_temp  = intersectCoord_final;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Add the neighboring matrix cells to the list for intersection check if it is the first try
                % but no intersection ocuurs, so another one must be checked
                if ( ( isempty(intersectCoord_temp) ) || ( size(intersectCoord_temp,2) < 3 ) ) && Im_count == 2
                    Im_next = Index_Matrix_3D( NX,NY,NZ,max(i-1,1 ),j          ,k           );  if ~ismember(Im_next,Im_List),  Im_List = [Im_List;Im_next];  end
                    Im_next = Index_Matrix_3D( NX,NY,NZ,min(i+1,NX),j          ,k           );  if ~ismember(Im_next,Im_List),  Im_List = [Im_List;Im_next];  end
                    Im_next = Index_Matrix_3D( NX,NY,NZ,i          ,max(j-1,1 ),k           );  if ~ismember(Im_next,Im_List),  Im_List = [Im_List;Im_next];  end
                    Im_next = Index_Matrix_3D( NX,NY,NZ,i          ,min(j+1,NY),k           );  if ~ismember(Im_next,Im_List),  Im_List = [Im_List;Im_next];  end
                    Im_next = Index_Matrix_3D( NX,NY,NZ,i          ,j          ,max(k-1,1 ) );  if ~ismember(Im_next,Im_List),  Im_List = [Im_List;Im_next];  end
                    Im_next = Index_Matrix_3D( NX,NY,NZ,i          ,j          ,min(k+1,NZ) );  if ~ismember(Im_next,Im_List),  Im_List = [Im_List;Im_next];  end
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Assigning the index of matrix cell to the array
                if ( isempty(intersectCoord_temp) ),  continue;  end
                
                % Add the neighboring cells to the list for intersection check
                Im_next = Index_Matrix_3D( NX,NY,NZ,max(i-1,1 ),j          ,k           );  if ~ismember(Im_next,Im_List),  Im_List = [Im_List;Im_next];  end
                Im_next = Index_Matrix_3D( NX,NY,NZ,min(i+1,NX),j          ,k           );  if ~ismember(Im_next,Im_List),  Im_List = [Im_List;Im_next];  end
                Im_next = Index_Matrix_3D( NX,NY,NZ,i          ,max(j-1,1 ),k           );  if ~ismember(Im_next,Im_List),  Im_List = [Im_List;Im_next];  end
                Im_next = Index_Matrix_3D( NX,NY,NZ,i          ,min(j+1,NY),k           );  if ~ismember(Im_next,Im_List),  Im_List = [Im_List;Im_next];  end
                Im_next = Index_Matrix_3D( NX,NY,NZ,i          ,j          ,max(k-1,1 ) );  if ~ismember(Im_next,Im_List),  Im_List = [Im_List;Im_next];  end
                Im_next = Index_Matrix_3D( NX,NY,NZ,i          ,j          ,min(k+1,NZ) );  if ~ismember(Im_next,Im_List),  Im_List = [Im_List;Im_next];  end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Ordering the intersection Points (either clockwise or counter clockwise)
                if size(intersectCoord_temp,2) > 2
                    distanceToZero = sqrt( sum( intersectCoord_temp.^2 , 1 ) );

                    % Setting the 1st Point
                    [ ~ , index ] = min( distanceToZero );
                    intersectCoord_final      = zeros ( size(intersectCoord_temp,1) , size(intersectCoord_temp,2) );
                    intersectCoord_final(:,1) = intersectCoord_temp(:,index);
                    intersectCoord_temp(:,index) = [];

                    % Setting the 2nd Point
                    distanceToRef = sqrt( sum ( ( intersectCoord_temp - intersectCoord_final(:,1) ).^2 , 1 ) );
                    [ ~ , index ] = min( distanceToRef );
                    intersectCoord_final(:,2) = intersectCoord_temp(:,index);
                    intersectCoord_temp(:,index) = [];

                    % Setting the rest of the Points based on the maximum angle between two close segments
                    for m = 2 : length( distanceToZero ) - 2
                        arcCos_Theta = ones( size( intersectCoord_temp , 2 ) , 1 );
                        for n = 1 : size( intersectCoord_temp , 2 )
                            arcCos_Theta(n) = dot( intersectCoord_final(:,m) - intersectCoord_final(:,m-1) , ...
                                                   intersectCoord_final(:,m) - intersectCoord_temp (:,n  ) ) ...
                                              / norm( intersectCoord_final(:,m) - intersectCoord_final(:,m-1) ) ...
                                              / norm( intersectCoord_final(:,m) - intersectCoord_temp (:,n  ) );
                        end
                        [ ~ , index ] = min( arcCos_Theta );
                        intersectCoord_final(:,m+1)  = intersectCoord_temp(:,index);
                        intersectCoord_temp(:,index) = [];
                    end
                    intersectCoord_final(:,end) = intersectCoord_temp;
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Writing the final intersection Points into the fracture data structure
                if ( size(intersectCoord_final,2) < 3 ),  continue;  end
                
                Index_matIntersect = Index_matIntersect + 1;
                Frac(f).intersectCoord_matCell{If}{Index_matIntersect,1} = Im; 
                Frac(f).      areaFrac_matCell{If}{Index_matIntersect,1} = Im;
                Frac(f).       aveDist_matCell{If}{Index_matIntersect,1} = Im;
                Frac(f).         T_Geo_matCell{If}{Index_matIntersect,1} = Im;
                
                Frac(f).intersectCoord_matCell{If}{Index_matIntersect,2} = intersectCoord_final;

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Calculating the area fraction
                Frac(f).areaFrac_matCell{If}{Index_matIntersect,2} = 0;
%                 if ( size(intersectCoord_final,2) < 3 )
%                     Frac(f).aveDist_matCell{If}{Index_matIntersect,2} = 1;  % If areaFrac_matCell is zero, make aveDist_matcube 1;
%                     Frac(f).  T_Geo_matCell{If}{Index_matIntersect,2} = 0;  % If areaFrac_matCell is zero, make   T_Geo_matcube 0;
%                     continue;
%                 end
                
                for n = 2 : size( intersectCoord_final , 2 )  - 1
                    Frac(f).areaFrac_matCell{If}{Index_matIntersect,2} = Frac(f).areaFrac_matCell{If}{Index_matIntersect,2} + ...
                        Triangle_Area_3D( intersectCoord_final(:,1) , intersectCoord_final(:,n) , intersectCoord_final(:,n+1) );
                end
                
                % Doubling the area fraction if the cell is not coplanar to any matrix cube edges
                if areCoplanar == 0
                    Frac(f).areaFrac_matCell{If}{Index_matIntersect,2} = Frac(f).areaFrac_matCell{If}{Index_matIntersect,2} * 2;
                end
                Frac(f).areaFrac_matCell_sum = Frac(f).areaFrac_matCell_sum + Frac(f).areaFrac_matCell{If}{Index_matIntersect,2};

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Obtaining the average distance between the matrix cube cell and the fracture plate cell
                Frac(f).aveDist_matCell{If}{Index_matIntersect,2} = 0; 
                N_cube = 06;
                Xim_cube = linspace( Xim(i), Xim(i+1), N_cube );
                Yim_cube = linspace( Yim(j), Yim(j+1), N_cube );
                Zim_cube = linspace( Zim(k), Zim(k+1), N_cube );
                for kk = 1 : N_cube
                    for jj = 1 : N_cube
                        Frac(f).aveDist_matCell{If}{Index_matIntersect,2} = Frac(f).aveDist_matCell{If}{Index_matIntersect,2} + ...
                            sum ( abs( Frac(f).Equation.a * Xim_cube + Frac(f).Equation.b * Yim_cube(jj) + Frac(f).Equation.c * Zim_cube(kk) - Frac(f).Equation.d ) );
                    end
                end
                Frac(f).aveDist_matCell{If}{Index_matIntersect,2} = Frac(f).aveDist_matCell{If}{Index_matIntersect,2} ...
                                                                  / sqrt( Frac(f).Equation.a^2 + Frac(f).Equation.b^2 + Frac(f).Equation.c^2 ) / N_cube^3 ;
                                                              
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Obtaining the geometrical transmissibility between the matrix cube cell and the fracture plate cell ( CI = Af/<d> )
                a = Frac(f).areaFrac_matCell{If}{Index_matIntersect,2} / Frac(f).aveDist_matCell{If}{Index_matIntersect,2};
                Frac(f).T_Geo_matCell{If}{Index_matIntersect,2} = a;
                                                  
            end % End of while-loop
            
            % Obtaining the overlapped cells due to the intersection line between fractures f and g
            max_dist = sqrt( Frac(f).D_Length_AB^2 + Frac(f).D_Width_AD^2 )/2;
            for g = 1 : length(Frac)
                if ~isempty(Frac(f).intersectCoord_fracObj{g})
                    p1 = Frac(f).intersectCoord_fracObj{g}(:,1);
                    p2 = Frac(f).intersectCoord_fracObj{g}(:,2);
                    Frac(f).Overlap_frac{g}((j_f-1)*Frac(f).N_Length_AB+i_f) = ( norm( cross(Frac(f).CellCenterCoordsV2((j_f-1)*Frac(f).N_Length_AB+i_f,:)'-p1 , p2-p1) ) / norm(p2-p1) ) <= max_dist;
                end
            end

        end % End of fractre 1st for-loop
    end % End of fractre 2nd for-loop
    
end % End of main fracture for-loop
fprintf('---------------------------------------------------------\n');
%load 'D:\SURFdrive\Simulation\DARSim2\Input_Desktop_Dell\ImmHomo\FracGen_3D_Workspace_Temp_All.mat';
%% Intersections of fracture cells of distinct fractures with each other (if any)
fprintf('Obtaining fracture - fracture connectivities:\n');
almostZero = almostZero * 10;
for f = 1 : length(Frac)
    for g = f+1 : length(Frac)
        % Finding the intersections between the cells of distinct fractures
        if ~isempty(Frac(f).intersectCoord_fracObj{g})
            fprintf('Fracture %02d <---> %02d\n',f,g);
            Frac(f).intersectCoord_fracCell{g} = cell( Frac(f).N_Length_AB * Frac(f).N_Width_AD , 1 );
            Frac(g).intersectCoord_fracCell{f} = cell( Frac(g).N_Length_AB * Frac(g).N_Width_AD , 1 );

            Frac(f).      areaFrac_fracCell{g} = cell( Frac(f).N_Length_AB * Frac(f).N_Width_AD , 1 );
            Frac(f).       aveDist_fracCell{g} = cell( Frac(f).N_Length_AB * Frac(f).N_Width_AD , 1 );
            Frac(f).  areaFrac_fracCell_sum{g} = 0;
             
            Frac(g).      areaFrac_fracCell{f} = cell( Frac(g).N_Length_AB * Frac(g).N_Width_AD , 1 );
            Frac(g).       aveDist_fracCell{f} = cell( Frac(g).N_Length_AB * Frac(g).N_Width_AD , 1 );     
            Frac(g).  areaFrac_fracCell_sum{f} = 0;
            
            % Loop for fracture f
            If_List = find(Frac(f).Overlap_frac{g});
            for nf = 1:length(If_List)
                If = If_List(nf);
                i_f = mod( If                           , Frac(f).N_Length_AB )   ;   if ( i_f==0 ),   i_f = Frac(f).N_Length_AB;   end
                j_f = mod( (If-i_f)/Frac(f).N_Length_AB , Frac(f).N_Width_AD  ) +1;   if ( j_f==0 ),   j_f = Frac(f).N_Width_AD ;   end
                if Index_Fracture_2D( Frac(f).N_Length_AB , Frac(f).N_Width_AD , i_f , j_f ) ~= If
                    error('i_f,j_f are not correspondent with If. Check the formula again!');
                end
                    
                % Corner Points of the 1st fracture cell
                Plane_fracCell_1.PointA = [ Frac(f).GridCoords(i_f  ,j_f  ,1) ; Frac(f).GridCoords(i_f  ,j_f  ,2) ; Frac(f).GridCoords(i_f  ,j_f  ,3) ];
                Plane_fracCell_1.PointB = [ Frac(f).GridCoords(i_f  ,j_f+1,1) ; Frac(f).GridCoords(i_f  ,j_f+1,2) ; Frac(f).GridCoords(i_f  ,j_f+1,3) ];
                Plane_fracCell_1.PointC = [ Frac(f).GridCoords(i_f+1,j_f+1,1) ; Frac(f).GridCoords(i_f+1,j_f+1,2) ; Frac(f).GridCoords(i_f+1,j_f+1,3) ];
                Plane_fracCell_1.PointD = [ Frac(f).GridCoords(i_f+1,j_f  ,1) ; Frac(f).GridCoords(i_f+1,j_f  ,2) ; Frac(f).GridCoords(i_f+1,j_f  ,3) ];
                
                % Loop for fracture g
                Index_frac_Intersect  = 0;
                Index_frac_Af_aveDist = 0;
                Ig_List = find(Frac(g).Overlap_frac{f});
                for ng = 1:length(Ig_List)
                    Ig = Ig_List(ng);
                    i_g = mod( Ig                           , Frac(g).N_Length_AB )   ;   if ( i_g==0 ),   i_g = Frac(g).N_Length_AB;   end
                    j_g = mod( (Ig-i_g)/Frac(g).N_Length_AB , Frac(g).N_Width_AD  ) +1;   if ( j_g==0 ),   j_g = Frac(g).N_Width_AD ;   end
                    if Index_Fracture_2D( Frac(g).N_Length_AB , Frac(g).N_Width_AD , i_g , j_g ) ~= Ig
                        error('i_g,j_g are not correspondent with Ig. Check the formula again!');
                    end
                    
                    % Corner Points of the 2nd fracture cell
                    Plane_fracCell_2.PointA = [ Frac(g).GridCoords(i_g  ,j_g  ,1) ; Frac(g).GridCoords(i_g  ,j_g  ,2) ; Frac(g).GridCoords(i_g  ,j_g  ,3) ];
                    Plane_fracCell_2.PointB = [ Frac(g).GridCoords(i_g  ,j_g+1,1) ; Frac(g).GridCoords(i_g  ,j_g+1,2) ; Frac(g).GridCoords(i_g  ,j_g+1,3) ];
                    Plane_fracCell_2.PointC = [ Frac(g).GridCoords(i_g+1,j_g+1,1) ; Frac(g).GridCoords(i_g+1,j_g+1,2) ; Frac(g).GridCoords(i_g+1,j_g+1,3) ];
                    Plane_fracCell_2.PointD = [ Frac(g).GridCoords(i_g+1,j_g  ,1) ; Frac(g).GridCoords(i_g+1,j_g  ,2) ; Frac(g).GridCoords(i_g+1,j_g  ,3) ];
                    
                    [~, intersectPoints] = Plane_Seg_Intersect_3D( Plane_fracCell_1 , Plane_fracCell_2 , almostZero );
                    
                    intersectCoord_fracCell_final = [];
                    
                    if ~isempty(intersectPoints)
                        for nr = 1 : size(intersectPoints,2)
                            
                            % Removing the Points that are not inside the fracture rectangle
                            isInside1 = Is_Point_Inside_Rectangle_3D( Plane_fracCell_1 , intersectPoints(:,nr) , almostZero );
                            isInside2 = Is_Point_Inside_Rectangle_3D( Plane_fracCell_2 , intersectPoints(:,nr) , almostZero );
                            if ( isInside1 == 0 ) || ( isInside2 == 0 )
                                intersectPoints(:,nr) = [ dummy ; dummy ; dummy ];
                                continue;
                            end
                            
                            % Removing the Points that are repeated
                            subtracted = ( intersectPoints - intersectPoints(:,nr) );
                            subtracted = sqrt( subtracted(1,:).^2 + subtracted(2,:).^2 + subtracted(3,:).^2 );
                            if ~isempty( find( subtracted([1:nr-1 , nr+1:end]) < almostZero , 1 ) )
                                intersectPoints(:,nr) = [ dummy ; dummy ; dummy ];
                                continue;
                            end
                            
                            if intersectPoints(1,nr) ~= dummy
                                intersectCoord_fracCell_final = [ intersectCoord_fracCell_final , intersectPoints(:,nr) ];
                            end
                        end
                        
                        if ( size( intersectCoord_fracCell_final , 2 ) > 2 )
                            error('More than two intersection Points exist between two intersecting fracture cells!');
                        end
                        
                        % Writing intersection coordinates data into the structure
                        if ~isempty( intersectCoord_fracCell_final )
                            Index_frac_Intersect = Index_frac_Intersect + 1;
                            Frac(f).intersectCoord_fracCell{g}{If}{Index_frac_Intersect,1} = Ig;
                            Frac(f).intersectCoord_fracCell{g}{If}{Index_frac_Intersect,2} = intersectCoord_fracCell_final;
                            Frac(g).intersectCoord_fracCell{f}{Ig}{Index_frac_Intersect,1} = If;
                            Frac(g).intersectCoord_fracCell{f}{Ig}{Index_frac_Intersect,2} = intersectCoord_fracCell_final;
                        end
                        
                        % Obtaining the area fractions and average distances and writing them into structure
                        if ( size( intersectCoord_fracCell_final , 2 ) == 2 )
                            Frac(f).NumOf_fracCellConn(If) = Frac(f).NumOf_fracCellConn(If) + 1;
                            Frac(g).NumOf_fracCellConn(Ig) = Frac(g).NumOf_fracCellConn(Ig) + 1;
                            Frac(f).areaFrac_fracCell{g}{If}{Frac(f).NumOf_fracCellConn(If),1} = Ig;
                            Frac(f). aveDist_fracCell{g}{If}{Frac(f).NumOf_fracCellConn(If),1} = Ig;
                            Frac(g).areaFrac_fracCell{f}{Ig}{Frac(g).NumOf_fracCellConn(Ig),1} = If;
                            Frac(g). aveDist_fracCell{f}{Ig}{Frac(g).NumOf_fracCellConn(Ig),1} = If;
                            
                            [ Frac(f).areaFrac_fracCell{g}{If}{Frac(f).NumOf_fracCellConn(If),2} , Frac(f).aveDist_fracCell{g}{If}{Frac(f).NumOf_fracCellConn(If),2} , Collinearity1 ] = ...
                                Line_Plane_Connectivity_3D( Plane_fracCell_1 , intersectCoord_fracCell_final(:,1) , intersectCoord_fracCell_final(:,2) , almostZero );
                            [ Frac(g).areaFrac_fracCell{f}{Ig}{Frac(g).NumOf_fracCellConn(Ig),2} , Frac(g).aveDist_fracCell{f}{Ig}{Frac(g).NumOf_fracCellConn(Ig),2} , Collinearity2 ] = ...
                                Line_Plane_Connectivity_3D( Plane_fracCell_2 , intersectCoord_fracCell_final(:,1) , intersectCoord_fracCell_final(:,2) , almostZero );
                            
                            %                                     if Frac(f).areaFrac_fracCell{g}{If}{Frac(f).NumOf_fracCellConn(If),2} ~= Frac(g).areaFrac_fracCell{f}{Ig}{Frac(g).NumOf_fracCellConn(Ig),2}
                            %                                         error('WTF!');
                            %                                     end
                            
                            if ( Collinearity1 == 1 ),  Frac(g).areaFrac_fracCell{f}{Ig}{Frac(g).NumOf_fracCellConn(Ig),2} = Frac(g).areaFrac_fracCell{f}{Ig}{Frac(g).NumOf_fracCellConn(Ig),2} /2;  end
                            if ( Collinearity2 == 1 ),  Frac(f).areaFrac_fracCell{g}{If}{Frac(f).NumOf_fracCellConn(If),2} = Frac(f).areaFrac_fracCell{g}{If}{Frac(f).NumOf_fracCellConn(If),2} /2;  end
                            
                            T_Geo_1 = Frac(f).areaFrac_fracCell{g}{If}{Frac(f).NumOf_fracCellConn(If),2} / Frac(f).aveDist_fracCell{g}{If}{Frac(f).NumOf_fracCellConn(If),2};
                            T_Geo_2 = Frac(g).areaFrac_fracCell{f}{Ig}{Frac(g).NumOf_fracCellConn(Ig),2} / Frac(g).aveDist_fracCell{f}{Ig}{Frac(g).NumOf_fracCellConn(Ig),2};
                            T_Geo_Harmonic = 2*T_Geo_1*T_Geo_2 / (T_Geo_1 + T_Geo_2);
                            
                            Frac(f).T_Geo_fracCell{g}{If}{Frac(f).NumOf_fracCellConn(If),2} = T_Geo_Harmonic;
                            Frac(g).T_Geo_fracCell{f}{Ig}{Frac(g).NumOf_fracCellConn(Ig),2} = T_Geo_Harmonic;
                            
                            Frac(f).areaFrac_fracCell_sum{g} = Frac(f).areaFrac_fracCell_sum{g} + Frac(f).areaFrac_fracCell{g}{If}{Frac(f).NumOf_fracCellConn(If),2};
                            Frac(g).areaFrac_fracCell_sum{f} = Frac(g).areaFrac_fracCell_sum{f} + Frac(g).areaFrac_fracCell{f}{Ig}{Frac(g).NumOf_fracCellConn(Ig),2};
                            
                        end
                        
                        
                    end
                    
                end
                
            end
            
        end
 
    end
end

%% Obtaining the number frac-frac connectivities for each fracture cell
% for f = 1 : length(Frac)
%     Frac(f).NumOf_fracCellConn = zeros( Frac(f).N_Length_AB * Frac(f).N_Width_AD , 1 );
%     for g = 1 : length(Frac)  
%         if f ~= g
%             for If = 1 : length( Frac(f).areaFrac_fracCell{g} )
%                 if ~isempty( Frac(f).areaFrac_fracCell{g}{If} )
%                     temp = Frac(f).areaFrac_fracCell{g}{If};
%                     temp(all(cellfun(@isempty,temp),2), : ) = [];
%                     Frac(f).areaFrac_fracCell{g}{If} = temp;
%                 end
%                 Frac(f).NumOf_fracCellConn(If) = Frac(f).NumOf_fracCellConn(If) + size( Frac(f).areaFrac_fracCell{g}{If},1 );
%             end
%         end
%     end
% end
%% End of Function
end