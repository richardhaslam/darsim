% class fracturer for DARSim2FracGen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Mousa HosseiniMehr
%TU Delft
%Created: 2017-07-12
%Last modified: 2019-02-08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef fracture_FracGen < planeSegment_FracGen
    properties
        Aperture
        Porosity
        Permeability
        Area_Total
        n_vec_plot
        N_Length_AB
        N_Width_AD
        N_Total
        ADM
        D_Length_AB
        D_Width_AD
        D_Length_AB_vec
        D_Width_AD_vec
        CellCenterCoordsV1
        CellCenterCoordsV2
        GridCoords
        intersectCoord_rock
        Af_rock_EDFM
        Af_rock_EDFM_sum
        avgDist_rock_EDFM
        CI_rock_EDFM
        Af_rock_pEDFM
        avgDist_rock_pEDFM
        CI_rock_pEDFM
        MatrixOverlapList = [];
        intersectCoord_fracturePlate
        intersectCoord_frac
        Overlap_frac
        Af_frac_EDFM
        Af_frac_EDFM_sum
        avgDist_frac_EDFM
        CI_frac_EDFM
        NumOf_fracConn_EDFM
        Af_frac_pEDFM
        avgDist_frac_pEDFM
        CI_frac_pEDFM
        NumOf_fracConn_pEDFM
        alpha_Tx
        alpha_Ty
        % alpha is a the fraction of blockage (from 0 to 1) between
        % neighboring fracture cells due to existence of other intersecting
        % fracture planes when considering pEDFM method. The alpha=0 means 
        % no blockage and alpha=1 means total blockage.
    end
    methods
        function obj = fracture_FracGen(PointA,PointB,PointC,PointD)
            if nargin==4
                obj.InitializePlaneSegment(obj,PointA,PointB,PointC,PointD);
            end
        end
        function BuildFracture(obj, Frac_Input, ReservoirGrid, Epsilon, f)
            % Dimensions of Reservoir in each direction
            LX = ReservoirGrid.LX;
            LY = ReservoirGrid.LY;
            LZ = ReservoirGrid.LZ;
            
            if isnan( Frac_Input.CenterCoord(1) )
                % Fracture construction with three Points input type:
                % The three vectors connecting 3 corners of the fracture plate (forming a triangle)
                Vec12 = Frac_Input.CornerCoords(:,2) - Frac_Input.CornerCoords(:,1);
                Vec23 = Frac_Input.CornerCoords(:,3) - Frac_Input.CornerCoords(:,2);
                Vec31 = Frac_Input.CornerCoords(:,3) - Frac_Input.CornerCoords(:,1);
                
                % Finding which of the three Points is the corner with the 90 deg angle
                dotProduct = [ dot( Vec31 , Vec12 ) ; dot( Vec12 , Vec23 ) ; dot( Vec23 , Vec31 ) ];
                PointRightAngle   = find( abs(dotProduct) <= Epsilon , 1 );
                PointSalientAngle = find( abs(dotProduct) >  Epsilon     );
                if isempty(PointRightAngle),  error('Fracture # %02d is not a rectangle! Please check the input file.',f);  end
                
                % Setting attributes to the corners (A is the corner with 90 deg angle)
                PointA = Frac_Input.CornerCoords( : , PointRightAngle      );
                PointB = Frac_Input.CornerCoords( : , PointSalientAngle(1) );
                PointD = Frac_Input.CornerCoords( : , PointSalientAngle(2) );
                
                % The central Point of the plate
                PointM = ( PointB + PointD ) /2;
                
                % Obtaining the 4th corner (by using the centeral Point)
                PointC = PointA + 2* ( PointM - PointA );
            else
                
                % Fracture construction with central Point and angles input type
                PointM = Frac_Input.CenterCoord;
                
                % Obtaining the corners before rotation
                PointA_noRot       = [  PointM(1) - Frac_Input.Length /2  ;  PointM(2) + Frac_Input.Width /2  ;  PointM(3)  ];
                PointB_noRot       = [  PointM(1) + Frac_Input.Length /2  ;  PointM(2) + Frac_Input.Width /2  ;  PointM(3)  ];
                PointC_noRot       = [  PointM(1) + Frac_Input.Length /2  ;  PointM(2) - Frac_Input.Width /2  ;  PointM(3)  ];
                PointD_noRot       = [  PointM(1) - Frac_Input.Length /2  ;  PointM(2) - Frac_Input.Width /2  ;  PointM(3)  ];
                
                % Rotating around Z axis with horizontal rotation angle
                RotAxis_Hor        = [ 0 ; 0 ; 1 ];
                RotAngle_Hor_Rad   = Frac_Input.RotAngleAlongZ * pi / 180;
                PointA_Rot_Hor     = Rotate_Point_Around_Line_3D( PointA_noRot , RotAxis_Hor , PointM , RotAngle_Hor_Rad );
                PointB_Rot_Hor     = Rotate_Point_Around_Line_3D( PointB_noRot , RotAxis_Hor , PointM , RotAngle_Hor_Rad );
                PointC_Rot_Hor     = Rotate_Point_Around_Line_3D( PointC_noRot , RotAxis_Hor , PointM , RotAngle_Hor_Rad );
                PointD_Rot_Hor     = Rotate_Point_Around_Line_3D( PointD_noRot , RotAxis_Hor , PointM , RotAngle_Hor_Rad );
                
                % Rotating around "BMC - DMA" (along central length)
                PointBMC           = ( PointB_Rot_Hor + PointC_Rot_Hor ) /2;
                PointDMA           = ( PointD_Rot_Hor + PointA_Rot_Hor ) /2;
                RotAxis_Length     = PointDMA - PointBMC;
                RotAxis_Length_Rad  = Frac_Input.RotAngleAlongL * pi / 180;
                PointA_Rot_Length   = Rotate_Point_Around_Line_3D( PointA_Rot_Hor , RotAxis_Length , PointM , RotAxis_Length_Rad );
                PointB_Rot_Length   = Rotate_Point_Around_Line_3D( PointB_Rot_Hor , RotAxis_Length , PointM , RotAxis_Length_Rad );
                PointC_Rot_Length   = Rotate_Point_Around_Line_3D( PointC_Rot_Hor , RotAxis_Length , PointM , RotAxis_Length_Rad );
                PointD_Rot_Length   = Rotate_Point_Around_Line_3D( PointD_Rot_Hor , RotAxis_Length , PointM , RotAxis_Length_Rad );
                
                % Rotating around "AMB - CMD" (along central width)
                PointAMB           = ( PointA_Rot_Length + PointB_Rot_Length ) /2;
                PointCMD           = ( PointC_Rot_Length + PointD_Rot_Length ) /2;
                RotAxis_Width      = PointCMD - PointAMB;
                RotAxis_Width_Rad = Frac_Input.RotAngleAlongW * pi / 180;
                PointA     = Rotate_Point_Around_Line_3D( PointA_Rot_Length , RotAxis_Width , PointM , RotAxis_Width_Rad );
                PointB     = Rotate_Point_Around_Line_3D( PointB_Rot_Length , RotAxis_Width , PointM , RotAxis_Width_Rad );
                PointC     = Rotate_Point_Around_Line_3D( PointC_Rot_Length , RotAxis_Width , PointM , RotAxis_Width_Rad );
                PointD     = Rotate_Point_Around_Line_3D( PointD_Rot_Length , RotAxis_Width , PointM , RotAxis_Width_Rad );
            end
            
            %% Starting with discretization
            % Initializing the Plane Segment properties of the fracture
            obj.InitializePlaneSegment(PointA,PointB,PointC,PointD);
            
            % Checking if fracture fits inside the matrix box
            if ( min(obj.Points(1,:)) < 0 - Epsilon ) || ( max(obj.Points(1,:)) > LX + Epsilon ) || ...
               ( min(obj.Points(2,:)) < 0 - Epsilon ) || ( max(obj.Points(2,:)) > LY + Epsilon ) || ...
               ( min(obj.Points(3,:)) < 0 - Epsilon ) || ( max(obj.Points(3,:)) > LZ + Epsilon )
                error('Fracture #%02d does not fit in the matrix dimensions! Please check the input file.',f);
            end

            % The area of the fracture plane (each plane has two sides)
            obj.Area_Total = Triangle_Area_3D( obj.PointA , obj.PointB , obj.PointC ) * 2 * 2;
            
            % The normal vector of plate with its size normalized to 1/10 of the reservoir length (LX) (just for plotting)
            obj.n_vec_plot = obj.n_vec * LZ / 10;
            
            % Check if the fracture plate is paralel to any of X,Y,Z axes
            obj.CheckIfParallelToRefSurface(Epsilon);
            
            % Length (Length_AB) and width (Width_AD) of the fracture plate
            obj.Length_AB = norm( obj.AB_vec );
            obj.Width_AD  = norm( obj.AD_vec );
            
            % Adding ADM configuration
            obj.ADM = Frac_Input.ADM;
            
            % Discretizing the fracture plane
            obj.Discretize(Frac_Input, ReservoirGrid, f);
            
            % Aperture, Porosity and Permeability of The fracture plate
            obj.Aperture     = Frac_Input.Aperture;
            if isnan(obj.Aperture)
                obj.Aperture = mean([ReservoirGrid.DX,ReservoirGrid.DY])/10;
            end
            obj.Porosity     = Frac_Input.Porosity;
            obj.Permeability = Frac_Input.Permeability;
            
            %% Assigning Coordinates (x,y,z) of cell centers and cell interfaces of the fracture plane
            obj.AssignCellCoordinates;
        end
        function Discretize(obj, Frac_Input, ReservoirGrid, f)
            %%
            % Dimensions of Reservoir in each direction
            LX = ReservoirGrid.LX;
            LY = ReservoirGrid.LY;
            LZ = ReservoirGrid.LZ;
            % Number of reservoir grid cells in each direction
            NX = ReservoirGrid.NX;
            NY = ReservoirGrid.NY;
            NZ = ReservoirGrid.NZ;
            
            % Obtaining the number of grid cells at each direction
            if ~isnan( Frac_Input.GridResRatio )
                % Discretizing with "GridResRatio"
                obj.N_Length_AB = round( ( obj.Length_AB / max([LX,LY,LZ]) ) * max([NX,NY,NZ]) ^ Frac_Input.GridResRatio );
                obj.N_Width_AD  = round( ( obj.Width_AD  / max([LX,LY,LZ]) ) * max([NX,NY,NZ]) ^ Frac_Input.GridResRatio );
                if obj.ADM(1)
                    % There are two approaches for coarse grid constuctio:
                    % 1- Original, with identical primal coarse grids. In
                    % this case, fine grid number must be odd.
                    % 2- Coarse nodes on the corners and edges, resulting
                    % in unequal primal coarse grids (coarse grids on the
                    % corners and edges have smaller size as they are cut).
                    % In this case,  fine grid number must be even.
                    if mod(NX,2)==1 && mod(NY,2)==1
                        % Coarse grids will be constructed with original method
                        % Correcting N_Length_AB to match MMs level and coarsening ratio
                        if ( mod( obj.N_Length_AB , obj.ADM(3)^obj.ADM(2) )~=0 )
                            obj.N_Length_AB = obj.N_Length_AB - mod( obj.N_Length_AB , obj.ADM(3)^obj.ADM(2) );
                        end
                        % Correcting N_Width_AD to match MMs level and coarsening ratio
                        if ( mod( obj.N_Width_AD  , obj.ADM(4)^obj.ADM(2) )~=0 )
                            obj.N_Width_AD  = obj.N_Width_AD  - mod( obj.N_Width_AD  , obj.ADM(4)^obj.ADM(2) );
                        end
                        % At least 3 coarse grid cells at the highest coarsening level in each direction
                        obj.N_Length_AB = max( obj.N_Length_AB , 3*obj.ADM(3)^obj.ADM(2) );
                        obj.N_Width_AD  = max( obj.N_Width_AD  , 3*obj.ADM(4)^obj.ADM(2) );
                    else
                        % Coarse nodes will be put on the corners as well.
                        % Correcting N_Length_AB to match MMs level and coarsening ratio
                        if ( mod( obj.N_Length_AB-1 , obj.ADM(3)^obj.ADM(2) )~=0 )
                            obj.N_Length_AB = obj.N_Length_AB - mod( obj.N_Length_AB-1 , obj.ADM(3)^obj.ADM(2) );
                        end
                        % Correcting N_Width_AD to match MMs level and coarsening ratio
                        if ( mod( obj.N_Width_AD-1  , obj.ADM(4)^obj.ADM(2) )~=0 )
                            obj.N_Width_AD  = obj.N_Width_AD  - mod( obj.N_Width_AD  , obj.ADM(4)^obj.ADM(2) );
                        end
                        % At least 4 coarse grid cells at the highest coarsening level in each direction
                        obj.N_Length_AB = max( obj.N_Length_AB , 3*obj.ADM(3)^obj.ADM(2)+1 );
                        obj.N_Width_AD  = max( obj.N_Width_AD  , 3*obj.ADM(4)^obj.ADM(2)+1 );
                    end
                end

                % At least 3 grid cells in each direction
                obj.N_Length_AB = max( obj.N_Length_AB , 3);
                obj.N_Width_AD  = max( obj.N_Width_AD  , 3);
                
                % If reservoir has only 1 grid cell in Z direction 
                % (i.e., 2D reservoir), fracture must also have only
                % 1 grid cell along its width (i.e., 1D fracture)
                if NZ == 1,  obj.N_Width_AD = 1;  end
                
            else
                % Discretizing with "GridNumAlongL" and "GridNumAlongW"
                obj.N_Length_AB = Frac_Input.GridNumAlongL;
                obj.N_Width_AD  = Frac_Input.GridNumAlongW;
                
                if obj.ADM(1)
                    if mod(NX,2)==1 && mod(NY,2)==1
                        if mod( obj.N_Length_AB , obj.ADM(3)^obj.ADM(2) ) ~= 0
                            error('In fracture #%02d, the number of grid cells along length is not acceptable to the ADM properties!\nPlease check the input file.',f);
                        end
                        if mod( obj.N_Width_AD  , obj.ADM(4)^obj.ADM(2) ) ~= 0
                            error('In fracture #%02d, the number of grid cells along width is not acceptable to the ADM properties!\nPlease check the input file.',f);
                        end
                    else
                        if mod( obj.N_Length_AB-1 , obj.ADM(3)^obj.ADM(2) ) ~= 0
                            error('In fracture #%02d, the number of grid cells along length is not acceptable to the ADM properties!\nPlease check the input file.',f);
                        end
                        if mod( obj.N_Width_AD-1  , obj.ADM(4)^obj.ADM(2) ) ~= 0
                            error('In fracture #%02d, the number of grid cells along width is not acceptable to the ADM properties!\nPlease check the input file.',f);
                        end
                    end
                end
            end
            
            % Total number of grid cells in the fracture
            obj.N_Total = obj.N_Length_AB * obj.N_Width_AD;
            
            % The size of each fracture grid cell in AB and AD directions with their vectors
            obj.D_Length_AB     = obj.Length_AB / obj.N_Length_AB;
            obj.D_Width_AD      = obj.Width_AD  / obj.N_Width_AD;
            obj.D_Length_AB_vec = obj.AB_vec    / obj.N_Length_AB;
            obj.D_Width_AD_vec  = obj.AD_vec    / obj.N_Width_AD;
            
            % alpha is a the fraction of blockage (from 0 to 1) between
            % neighboring fracture cells due to existence of other intersecting
            % fracture planes when considering pEDFM method. The alpha=0 means
            % no blockage and alpha=1 means total blockage.
            obj.alpha_Tx = zeros( obj.N_Length_AB+1, obj.N_Width_AD   );
            obj.alpha_Ty = zeros( obj.N_Length_AB  , obj.N_Width_AD+1 );
        end
        function AssignCellCoordinates(obj)
            %% Coordinates (x,y,z) of cell centers and cell interfaces of the fracture plate
            obj.CellCenterCoordsV1  = zeros( obj.N_Length_AB   , obj.N_Width_AD   , 3 );
            obj.CellCenterCoordsV2  = zeros( obj.N_Length_AB * obj.N_Width_AD , 3 );
            obj.GridCoords = zeros( obj.N_Length_AB+1 , obj.N_Width_AD+1 , 3 );
            
            %% Coordinates of fracture cell centers
            for i = 1 : obj.N_Length_AB
                obj.CellCenterCoordsV1(i,:,1) = linspace( obj.PointA(1) + (i-1/2)*obj.D_Length_AB_vec(1) + obj.D_Width_AD_vec(1)/2 , ...
                    obj.PointD(1) + (i-1/2)*obj.D_Length_AB_vec(1) - obj.D_Width_AD_vec(1)/2 , obj.N_Width_AD );
                
                obj.CellCenterCoordsV1(i,:,2) = linspace( obj.PointA(2) + (i-1/2)*obj.D_Length_AB_vec(2) + obj.D_Width_AD_vec(2)/2 , ...
                    obj.PointD(2) + (i-1/2)*obj.D_Length_AB_vec(2) - obj.D_Width_AD_vec(2)/2 , obj.N_Width_AD );
                
                obj.CellCenterCoordsV1(i,:,3) = linspace( obj.PointA(3) + (i-1/2)*obj.D_Length_AB_vec(3) + obj.D_Width_AD_vec(3)/2 , ...
                    obj.PointD(3) + (i-1/2)*obj.D_Length_AB_vec(3) - obj.D_Width_AD_vec(3)/2 , obj.N_Width_AD );
            end
            
            %% Coordinates of fracture cell centers in vector form
            for j = 1 : obj.N_Width_AD
                obj.CellCenterCoordsV2( (j-1)*obj.N_Length_AB+1 : j*obj.N_Length_AB , : ) = obj.CellCenterCoordsV1(:,j,:);
            end
            
            %% Coordinates of fracture cell interfaces
            for i = 1 : obj.N_Length_AB+1
                obj.GridCoords(i,:,1) = linspace( obj.PointA(1) + (i-1)*obj.D_Length_AB_vec(1) , ...
                    obj.PointD(1) + (i-1)*obj.D_Length_AB_vec(1) , obj.N_Width_AD+1 );
                
                obj.GridCoords(i,:,2) = linspace( obj.PointA(2) + (i-1)*obj.D_Length_AB_vec(2) , ...
                    obj.PointD(2) + (i-1)*obj.D_Length_AB_vec(2) , obj.N_Width_AD+1 );
                
                obj.GridCoords(i,:,3) = linspace( obj.PointA(3) + (i-1)*obj.D_Length_AB_vec(3) , ...
                    obj.PointD(3) + (i-1)*obj.D_Length_AB_vec(3) , obj.N_Width_AD+1 );
            end
        end
        function PrintInfo(obj,f)
            %%
           
            fprintf('Fracture %2d: Dimension= %5.2f x %5.2f [m2] , Grid= %3.0f x %3.0f = %4.0f , ADM lvl= %1.0f\n', ...
                    f, obj.Length_AB, obj.Width_AD, obj.N_Length_AB, obj.N_Width_AD, obj.N_Length_AB*obj.N_Width_AD, obj.ADM(1)*obj.ADM(2) );
        end
    end
end