%  Prolongation builder MS pressure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 4 August 2017
%Last modified: 19 September 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef prolongation_builder_MSPressure < prolongation_builder
    properties
        C % correction functions
        BFUpdater
        Dimensions
        ADMmap
    end
    methods
        function obj = prolongation_builder_MSPressure(n, cf)
            obj@prolongation_builder(n)
            obj.R = cell(1, n);
            obj.P = cell(1, n);
            if cf(1,3) == 1 && cf(1,2) == 1
                obj.Dimensions = 1;
            elseif cf(1,3) == 1
                obj.Dimensions = 2;
            else
                obj.Dimensions = 3;
            end
            obj.ADMmap = adm_map(prod(cf, 2));
        end
        function BuildStaticOperators(obj, ProductionSystem, FluidModel, FineGrid, CrossConnections, maxLevel, CoarseGrid)
            % Initialise
            obj.R = cell(maxLevel(1), 1);
            obj.P = cell(maxLevel(1), 1);
            % Build Pressure system
            obj.BFUpdater.ConstructPressureSystem(ProductionSystem, FluidModel, FineGrid, CrossConnections);
            % Printing BF type on console
            if isprop(obj.BFUpdater,'BFtype')
                fprintf(char(strcat({'The Basis Functions are '},obj.BFUpdater.BFtype,'.\n')));
            end
            %Build static restriction operator (FV)
            disp('Building Restriction 1');
            obj.R{1} = obj.MsRestriction(FineGrid, CoarseGrid(:,1));
            % Build Prolongation operator
            disp('Building Prolongation 1');
            [obj.P{1}, obj.C{1}] = obj.BFUpdater.MsProlongation(FineGrid, CoarseGrid(:,1), obj.Dimensions);
            % Build tpfa coarse system of level 1 (with MsFE)
            obj.BFUpdater.UpdatePressureMatrix(obj.P{1}, CoarseGrid(:, 1));
            for x = 2:maxLevel(1)
                % Build static restriction operator (FV)
                disp(['Building Restriction ', num2str(x)]);
                obj.R{x} = obj.MsRestriction(CoarseGrid(:, x-1), CoarseGrid(:, x));
                % Build Prolongation operator
                disp(['Building Prolongation ', num2str(x)]);
                [obj.P{x}, obj.C{x}] = obj.BFUpdater.MsProlongation(CoarseGrid(:, x-1), CoarseGrid(:, x), obj.Dimensions);
                %Build tpfa coarse system of level x (with MsFE)
                obj.BFUpdater.UpdatePressureMatrix(obj.P{x}, CoarseGrid(:, x));
            end
        end
        function ADMProlp = ADMProlongation(obj, ADMGrid, GlobalGrids, ADMRest)
            % Pressure prolongation
            ADMProlp = 1;
            
            % Loop over the levels
            for level = ADMGrid.MaxLevel:-1:2
                Prolp = obj.LevelProlongation(ADMGrid, GlobalGrids(level), level);
                % Multiply by previous objects
                ADMProlp = Prolp * ADMProlp;
            end
            
            % Last prolongation is different coz I use fine-scale ordering
            Prolp = obj.LastProlongation(ADMGrid, GlobalGrids(1));
            
            % Multiply by previous objects
            ADMProlp = Prolp * ADMProlp;
        end
        function Prolp = LevelProlongation(obj, ADMGrid, FineGrid, level)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % For a given level the prolongation operator looks like this
            %       Nf     Nc
            %    ----       ----
            %    |      |      |
            % Nf |  I   |  0   |
            %    |______|______|
            %    |      |      |
            % Nx |      |      |
            %    |      |      |
            %    ----       ----
            
            
            % Update map for the next level
            obj.ADMmap.Update(ADMGrid, FineGrid, level);
            
            % 1. Build the object
            Prolp = sparse(sum(obj.ADMmap.Nf) + sum(obj.ADMmap.Nx), sum(obj.ADMmap.Nf) + sum(obj.ADMmap.Nc));
            
            % 2. Fill in top left
            Prolp(1:sum(obj.ADMmap.Nf), 1:sum(obj.ADMmap.Nf)) = speye(sum(obj.ADMmap.Nf));
            
            % 3. Fill in Bottom left
            Prolp(sum(obj.ADMmap.Nf) + 1 : end,  obj.ADMmap.Verteces) = obj.P{level}(obj.ADMmap.OriginalIndexNx, obj.ADMmap.OriginalIndexVerteces);
            
            % 4. Fill in Bottom right
            Prolp(sum(obj.ADMmap.Nf) + 1 :end, sum(obj.ADMmap.Nf) + 1 : end) = obj.P{level}(obj.ADMmap.OriginalIndexNx, obj.ADMmap.OriginalIndexNc);
        end
        function Prolp = LastProlongation(obj, ADMGrid, FineGrid)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % For the first level the prolongation operator looks like this
            %       Nf     Nl1
            %    ----       ----
            %    |      |      |
            %    |      |   0  |
            % Nl0|      |______|
            %    |      |      |
            %    |      |      |
            %    |      |      |
            %    ----       ----
            
            % Update map
            obj.ADMmap.Update(ADMGrid, FineGrid, 1);
            
            % 1. Build object
            Prolp = sparse(FineGrid.N, sum(obj.ADMmap.Nf) + sum(obj.ADMmap.Nc));
            
            % 2. Fill in FS verteces of level 1
            Prolp(:,  obj.ADMmap.Verteces) = obj.P{1}(:, obj.ADMmap.OriginalIndexVerteces);
            
            % 3. Fill in coarse-scale nodes
            Prolp(:, sum(obj.ADMmap.Nf) + 1 : end) = obj.P{1}(:, obj.ADMmap.OriginalIndexNc);
            
            % 4. Fill in fine-scale nodes first
            rows = obj.ADMmap.OriginalIndexNf';
            columns = 1:sum(obj.ADMmap.Nf);
            Prolp(rows, :) = 0; % if it s fine-scale already I get rid of useless fillings
            Prolp(sub2ind(size(Prolp), rows, columns)) = 1;
        end
        function UpdateProlongationOperator(obj, FineGrid, CoarseGrid, ProductionSystem)
            % for now I do not update pressure basis functions
        end
        function deltaP_w = StaticMultilevelPressureGuess(obj, ProductionSystem, FluidModel, FineGrid, CoarseGrid, CrossConnections)
            function Im = Index_Matrix( Nx,Ny,Nz, i,j,k )
                if (i>Nx),   error('Current "i" index exceeds Nx.');  end
                if (j>Ny),   error('Current "j" index exceeds Ny.');  end
                if (k>Nz),   error('Current "k" index exceeds Nz.');  end
                Im = (k-1)*Nx*Ny + (j-1)*Nx + i;
            end
            %% Initializing Variables
            obj.BFUpdater.MaxContrast = obj.BFUpdater.MaxContrast/obj.BFUpdater.MaxContrast;
            obj.BFUpdater.ConstructPressureSystem(ProductionSystem, FluidModel, FineGrid, CrossConnections);
            A = obj.BFUpdater.A;
            CF = [CoarseGrid(1).CoarseFactor(1);CoarseGrid(1).CoarseFactor(2);CoarseGrid(1).CoarseFactor(3)];
            Nmx = FineGrid(1).Nx;
            Nmy = FineGrid(1).Ny;
            Nmz = FineGrid(1).Nz;
            Nmt = Nmx*Nmy*Nmz;
            deltaP_w = zeros(size(A,1) , 1);
            P_init = ProductionSystem.Reservoir.State_old.Properties('P_1').Value;
            CIx = CoarseGrid(1).I(:,2);
            CIy = CoarseGrid(1).J(:,2);
            CIz = CoarseGrid(1).K(:,2);
            Inj = ProductionSystem.Wells.Inj;
            Prod = ProductionSystem.Wells.Prod;
            K = ProductionSystem.Reservoir.K;
            Sm = ProductionSystem.Reservoir.State.Properties('S_1').Value;
            Mob = FluidModel.ComputePhaseMobilities(Sm);
            
            %% Looping over the wells
            for w = 1:length(Inj)+length(Prod)
                if w<=length(Inj),  Well = Inj(w);
                else,               Well = Prod(w-length(Inj));   end

                %% Obtaining the frame of the local region
                % bx
                if ~any( CIx == Well.Coord(1,1) )
                    if     Well.Coord(1,1) < min(CIx),    bx = 1;
                    elseif Well.Coord(1,1) > max(CIx),    bx = max(CIx);
                    else,                                 bx = max(CIx(CIx<Well.Coord(1,1)));
                    end
                else
                    if     Well.Coord(1,1) == min(CIx),   bx = 1;
                    elseif Well.Coord(1,1) == max(CIx),   bx = max(CIx) - CF(1);
                    else,                                 bx = Well.Coord(1,1) - CF(1);
                    end
                end
                % ex
                if ~any( CIx == Well.Coord(1,2) )
                    if     Well.Coord(1,2) < min(CIx),    ex = min(CIx);
                    elseif Well.Coord(1,2) > max(CIx),    ex = Nmx;
                    else,                                 ex = min(CIx(CIx>Well.Coord(1,2)));
                    end
                else
                    if     Well.Coord(1,2) == min(CIx),   ex = min(CIx) + CF(1);
                    elseif Well.Coord(1,2) == max(CIx),   ex = Nmx;
                    else,                                 ex = Well.Coord(1,2) + CF(1);
                    end
                end
                
                % by
                if ~any( CIy == Well.Coord(2,1) )
                    if     Well.Coord(2,1) < min(CIy),    by = 1;
                    elseif Well.Coord(2,1) > max(CIy),    by = max(CIy);
                    else,                                 by = max(CIy(CIy<Well.Coord(2,1)));
                    end
                else
                    if     Well.Coord(2,1) == min(CIy),   by = 1;
                    elseif Well.Coord(2,1) == max(CIy),   by = max(CIy) - CF(2);
                    else,                                 by = Well.Coord(2,1) - CF(2);
                    end
                end
                % ey
                if ~any( CIy == Well.Coord(2,2) )
                    if     Well.Coord(2,2) < min(CIy),    ey = min(CIy);
                    elseif Well.Coord(2,2) > max(CIy),    ey = Nmy;
                    else,                                 ey = min(CIy(CIy>Well.Coord(2,2)));
                    end
                else
                    if     Well.Coord(2,2) == min(CIy),   ey = min(CIy) + CF(2);
                    elseif Well.Coord(2,2) == max(CIy),   ey = Nmy;
                    else,                                 ey = Well.Coord(2,2) + CF(2);
                    end
                end
                
                % bz
                if ~any( CIz == Well.Coord(3,1) )
                    if     Well.Coord(3,1) < min(CIz),    bz = 1;
                    elseif Well.Coord(3,1) > max(CIz),    bz = max(CIz);
                    else,                                 bz = max(CIz(CIz<Well.Coord(3,1)));
                    end
                else
                    if     Well.Coord(3,1) == min(CIz),   bz = 1;
                    elseif Well.Coord(3,1) == max(CIz),   bz = max(CIz) - CF(3);
                    else,                                 bz = Well.Coord(3,1) - CF(3);
                    end
                end
                % ez
                if ~any( CIz == Well.Coord(3,2) )
                    if     Well.Coord(3,2) < min(CIz),    ez = min(CIz);
                    elseif Well.Coord(3,2) > max(CIz),    ez = Nmz;
                    else,                                 ez = min(CIz(CIz>Well.Coord(3,2)));
                    end
                else
                    if     Well.Coord(3,2) == min(CIz),   ez = min(CIz) + CF(3);
                    elseif Well.Coord(3,2) == max(CIz),   ez = Nmz;
                    else,                                 ez = Well.Coord(3,2) + CF(3);
                    end
                end
                if Nmx == 1,  ex = bx;  end
                if Nmy == 1,  ey = by;  end
                if Nmz == 1,  ez = bz;  end 
                
                NxL = ex-bx+1;
                NyL = ey-by+1;
                NzL = ez-bz+1;
                NtL = NxL*NyL*NzL;
                
                frac = [];
                for k = bz : ez
                    for j = by : ey
                        for i = bx : ex
                            Im = Index_Matrix( Nmx,Nmy,Nmz , i ,j ,k  );
                            frac = [ frac ; find(A(Im,Nmt+1:end)'~=0)+Nmt ];
                        end
                    end
                end
                frac = unique(frac);
                Nf = length(frac);
                
                %% Local 3D Initialization
                AL3D = zeros(NxL*NyL*NzL+Nf  );
                qL3D = zeros(NxL*NyL*NzL+Nf,1);
                iLcx3D = [];  jLcy3D = [];  kLcz3D = [];       % Local indexes of inner coarse nodes
                for kL = 1:NzL
                    k = bz+kL-1;
                    for jL = 1:NyL
                        j = by+jL-1;
                        for iL = 1:NxL
                            i = bx+iL-1;
                            Im    = Index_Matrix( Nmx,Nmy,Nmz , i ,j ,k  );
                            ImL3D = Index_Matrix( NxL,NyL,NzL , iL,jL,kL );
                            if ( iL > 1 )
                                I_Left   = Index_Matrix( Nmx,Nmy,Nmz , i -1,j ,k  );
                                I_Left_L = Index_Matrix( NxL,NyL,NzL , iL-1,jL,kL );
                                AL3D( ImL3D , I_Left_L ) = A   ( Im    , I_Left );
                                AL3D( ImL3D , ImL3D    ) = AL3D( ImL3D , ImL3D  ) - A( Im , I_Left );
                            end
                            if ( iL < NxL )
                                I_Right   = Index_Matrix( Nmx,Nmy,Nmz , i +1,j ,k  );
                                I_Right_L = Index_Matrix( NxL,NyL,NzL , iL+1,jL,kL );
                                AL3D( ImL3D , I_Right_L ) = A   ( Im    , I_Right );
                                AL3D( ImL3D , ImL3D     ) = AL3D( ImL3D , ImL3D   ) - A( Im , I_Right );
                            end
                            if ( jL > 1 )
                                I_Back   = Index_Matrix( Nmx,Nmy,Nmz , i ,j -1,k  );
                                I_Back_L = Index_Matrix( NxL,NyL,NzL , iL,jL-1,kL );
                                AL3D( ImL3D , I_Back_L ) = A   ( Im    , I_Back );
                                AL3D( ImL3D , ImL3D    ) = AL3D( ImL3D , ImL3D  ) - A( Im , I_Back );
                            end
                            if ( jL < NyL )
                                I_Front   = Index_Matrix( Nmx,Nmy,Nmz , i ,j +1,k  );
                                I_Front_L = Index_Matrix( NxL,NyL,NzL , iL,jL+1,kL );
                                AL3D( ImL3D , I_Front_L ) = A   ( Im    , I_Front );
                                AL3D( ImL3D , ImL3D     ) = AL3D( ImL3D , ImL3D   ) - A( Im , I_Front );
                            end
                            if ( kL > 1 )
                                I_Bottom   = Index_Matrix( Nmx,Nmy,Nmz , i ,j ,k -1  );
                                I_Bottom_L = Index_Matrix( NxL,NyL,NzL , iL,jL,kL-1 );
                                AL3D( ImL3D , I_Bottom_L ) = A   ( Im    , I_Bottom );
                                AL3D( ImL3D , ImL3D      ) = AL3D( ImL3D , ImL3D    ) - A( Im , I_Bottom );
                            end
                            if ( kL < NzL )
                                I_Top   = Index_Matrix( Nmx,Nmy,Nmz , i ,j ,k +1  );
                                I_Top_L = Index_Matrix( NxL,NyL,NzL , iL,jL,kL+1 );
                                AL3D( ImL3D , I_Top_L ) = A   ( Im    , I_Top );
                                AL3D( ImL3D , ImL3D   ) = AL3D( ImL3D , ImL3D ) - A( Im , I_Top );
                            end
                            
                            % Adding matrix-fracture connections
                            AL3D( ImL3D , NtL+1:end ) = A( Im , frac' );
                            AL3D( NtL+1:end , ImL3D ) = A( frac , Im  );
                            AL3D( ImL3D , ImL3D ) = AL3D( ImL3D , ImL3D ) - sum( A( Im , frac' ) );
                            AL3D( sub2ind(size(AL3D), NtL+1:NtL+Nf, NtL+1:NtL+Nf) ) = AL3D( sub2ind(size(AL3D), NtL+1:NtL+Nf, NtL+1:NtL+Nf) ) - A( Im , frac' );
                            
                            % Getting the local indexes of inner coarse nodes
                            if any(CIx == i) && any(CIy == j) && any(CIz == k) 
                                if ( ~any(iLcx3D==iL) ) && ( i~=bx ) && ( i~=ex )
                                    iLcx3D = [iLcx3D , iL];
                                end
                                if ( ~any(jLcy3D==jL) ) && ( j~=by ) && ( j~=ey )
                                    jLcy3D = [jLcy3D , jL];
                                end
                                if ( ~any(kLcz3D==kL) ) && ( k~=bz ) && ( k~=ez )
                                    kLcz3D = [kLcz3D , kL];
                                end
                                
                                % Setting the value of matrix coarse nodes to 0.0
                                if ~any(Well.Cells == Im)
                                    AL3D( ImL3D , :     ) = 0;
                                    AL3D( ImL3D , ImL3D ) = 1;
                                    qL3D( ImL3D         ) = 0;
                                end
                            end
                            
                            % Setting the value at the current well to 1.0 in AL3D
                            if any(Well.Cells == Im)
                                if w <= length(Inj)
                                    AL3D( ImL3D , ImL3D ) = AL3D( ImL3D , ImL3D ) + Well.PI * sum(Well.Mob, 2) * K(Im, 1);
                                    qL3D( ImL3D , 1     ) = Well.PI * sum(Well.Mob, 2) * K(Im, 1) * (Well.p-P_init(Im));
                                else
                                    AL3D( ImL3D , ImL3D ) = AL3D( ImL3D , ImL3D ) + Well.PI * sum(Mob(Im,:), 2) * K(Im, 1);
                                    qL3D( ImL3D , 1     ) = Well.PI * sum(Mob(Im,:), 2) * K(Im, 1) * (Well.p-P_init(Im));
                                end
                            end
                        end
                    end
                end
                
                % Adding fracture-fracture connections
                for f = 1 : Nf
                    AL3D( NtL+f , [NtL+1:NtL+f-1 , NtL+f+1:NtL+Nf] ) = A( frac(f) , [frac(1:f-1)' , frac(f+1:Nf)'] );
                    AL3D( [NtL+1:NtL+f-1 , NtL+f+1:NtL+Nf] , NtL+f ) = A( frac(f) , [frac(1:f-1)' , frac(f+1:Nf)'] )';
                    AL3D( NtL+f , NtL+f ) = AL3D( NtL+f , NtL+f ) - sum( A( frac(f) , [frac(1:f-1)' , frac(f+1:Nf)'] ) );
                end
                
                %% Localization in 2D-z faces
                for L2Dz = 1:length(kLcz3D)
                    AL2Dz = zeros( NxL*NyL     );
                    qL2Dz = zeros( NxL*NyL , 1 );
                    k = bz+kLcz3D(L2Dz)-1;
                    iLcx2D = [];  jLcy2D = [];
                    % AL2Dz
                    for jL = 1 : NyL
                        j = by+jL-1;
                        for iL = 1 : NxL
                            i = bx+iL-1;
                            Im    = Index_Matrix( Nmx,Nmy,Nmz , i ,j ,k );
                            ImL2D = Index_Matrix( NxL,NyL,1   , iL,jL,1 );
                            if ( iL > 1 )
                                I_Left   = Index_Matrix( Nmx,Nmy,Nmz , i -1,j ,k );
                                I_Left_L = Index_Matrix( NxL,NyL,1   , iL-1,jL,1 );
                                AL2Dz( ImL2D , I_Left_L ) = A    ( Im    , I_Left );
                                AL2Dz( ImL2D , ImL2D    ) = AL2Dz( ImL2D , ImL2D  ) - A( Im , I_Left );
                            end
                            if ( iL < NxL )
                                I_Right   = Index_Matrix( Nmx,Nmy,Nmz , i +1,j ,k );
                                I_Right_L = Index_Matrix( NxL,NyL,1   , iL+1,jL,1 );
                                AL2Dz( ImL2D , I_Right_L ) = A    ( Im    , I_Right );
                                AL2Dz( ImL2D , ImL2D     ) = AL2Dz( ImL2D , ImL2D   ) - A( Im , I_Right );
                            end
                            if ( jL > 1 )
                                I_Back   = Index_Matrix( Nmx,Nmy,Nmz , i ,j -1,k );
                                I_Back_L = Index_Matrix( NxL,NyL,1   , iL,jL-1,1 );
                                AL2Dz( ImL2D , I_Back_L ) = A    ( Im    , I_Back );
                                AL2Dz( ImL2D , ImL2D    ) = AL2Dz( ImL2D , ImL2D  ) - A( Im , I_Back );
                            end
                            if ( jL < NyL )
                                I_Front   = Index_Matrix( Nmx,Nmy,Nmz , i ,j +1,k );
                                I_Front_L = Index_Matrix( NxL,NyL,1   , iL,jL+1,1 );
                                AL2Dz( ImL2D , I_Front_L ) = A    ( Im    , I_Front );
                                AL2Dz( ImL2D , ImL2D     ) = AL2Dz( ImL2D , ImL2D   ) - A( Im , I_Front );
                            end
                            % Getting the local indexes of inner coarse nodes in 2D
                            if any(CIx == i) && any(CIy == j)
                                if ( i~=bx ) && ( i~=ex ) && ( ~any(iLcx2D==iL) )
                                    iLcx2D = [iLcx2D , iL];
                                end
                                if ( j~=by ) && ( j~=ey ) && ( ~any(jLcy2D==jL) )
                                    jLcy2D = [jLcy2D , jL];
                                end
                                % Setting the value of matrix coarse nodes to 0.0
                                if ~any(Well.Cells == Im)
                                    AL2Dz( ImL2D , :     ) = 0;
                                    AL2Dz( ImL2D , ImL2D ) = 1;
                                    AL2Dz( ImL2D         ) = 0;
                                end
                            end
                            % Setting the value at the current well to 1.0 in AL2Dz
                            if any(Well.Cells == Im)
                                if w <= length(Inj)
                                    AL2Dz( ImL2D , ImL2D ) = AL2Dz( ImL2D , ImL2D ) + Well.PI * sum(Well.Mob, 2) * K(Im, 1);
                                    qL2Dz( ImL2D , 1     ) = Well.PI * sum(Well.Mob, 2) * K(Im, 1) * (Well.p-P_init(Im));
                                else
                                    AL2Dz( ImL2D , ImL2D ) = AL2Dz( ImL2D , ImL2D ) + Well.PI * sum(Mob(Im,:), 2) * K(Im, 1);
                                    qL2Dz( ImL2D , 1     ) = Well.PI * sum(Mob(Im,:), 2) * K(Im, 1) * (Well.p-P_init(Im));
                                end
                            end
                        end
                    end
                    % Localization in 1D-x
                    for L1Dy = 1 : length(jLcy2D)
                        AL1Dx = zeros( NxL     );
                        qL1Dx = zeros( NxL , 1 );
                        j = by+jLcy2D(L1Dy)-1;
                        % AL1Dx
                        for iL = 1 : NxL
                            i = bx+iL-1;
                            Im = Index_Matrix( Nmx,Nmy,Nmz , i,j,k );
                            if ( iL > 1 )
                                I_Left = Index_Matrix( Nmx,Nmy,Nmz , i-1,j ,k );
                                AL1Dx( iL , iL-1 ) = A( Im , I_Left );
                                AL1Dx( iL , iL   ) = AL1Dx( iL , iL ) - A( Im , I_Left );
                            end
                            if ( iL < NxL )
                                I_Right = Index_Matrix( Nmx,Nmy,Nmz , i+1,j ,k );
                                AL1Dx( iL , iL+1 ) = A( Im , I_Right );
                                AL1Dx( iL , iL   ) = AL1Dx( iL , iL ) - A( Im , I_Right );
                            end
                            % Setting the value of matrix coarse nodes to 0.0
                            if ( any(CIx == i) )
                                if ~any(Well.Cells == Im)
                                    AL1Dx( iL , :   ) = 0;
                                    AL1Dx( iL , iL  ) = 1;
                                    qL1Dx( iL       ) = 0;
                                end
                            end
                            % Setting the value at the current well to 1.0 in AL1Dx
                            if any(Well.Cells == Im)
                                if w <= length(Inj)
                                    AL1Dx( iL , iL ) = AL1Dx( iL , iL ) + Well.PI * sum(Well.Mob, 2) * K(Im, 1);
                                    qL1Dx( iL , 1  ) = Well.PI * sum(Well.Mob, 2) * K(Im, 1) * (Well.p-P_init(Im));
                                else
                                    AL1Dx( iL , iL ) = AL1Dx( iL , iL ) + Well.PI * sum(Mob(Im,:), 2) * K(Im, 1);
                                    qL1Dx( iL , 1  ) = Well.PI * sum(Mob(Im,:), 2) * K(Im, 1) * (Well.p-P_init(Im));
                                end
                            end
                        end
                        % Dirichlet boundary condition at left
                        if bx > 1
                            AL1Dx( 1 , : ) = 0;
                            AL1Dx( 1 , 1 ) = 1;
                        end
                        % Dirichlet boundary condition at right
                        if ex < Nmx
                            AL1Dx( NxL , :   ) = 0;
                            AL1Dx( NxL , NxL ) = 1;
                        end
                        % Solving 1D localization on AL1Dx
                        pL1Dx = AL1Dx\qL1Dx;
                        % Assigning pL1Dx for Dirichlet Boundary Condition in AL2Dz
                        for iL = 1 : NxL
                            ImL2D = Index_Matrix( NxL,NyL,1, iL ,jLcy2D(L1Dy),1 );
                            AL2Dz( ImL2D , :     ) = 0;
                            AL2Dz( ImL2D , ImL2D ) = 1;
                            qL2Dz( ImL2D , 1     ) = pL1Dx(iL);
                        end
                    end
                    % Localization in 1D-y
                    for L1Dx = 1 : length(iLcx2D)
                        AL1Dy = zeros( NyL     );
                        qL1Dy = zeros( NyL , 1 );
                        i = bx+iLcx2D(L1Dx)-1;
                        % AL1Dy
                        for jL = 1 : NyL
                            j = by+jL-1;
                            Im = Index_Matrix( Nmx,Nmy,Nmz , i,j,k );
                            if ( jL > 1 )
                                I_Back = Index_Matrix( Nmx,Nmy,Nmz , i,j-1,k );
                                AL1Dy( jL , jL-1 ) = A( Im , I_Back );
                                AL1Dy( jL , jL   ) = AL1Dy( jL , jL ) - A( Im , I_Back );
                            end
                            if ( jL < NyL )
                                I_Front = Index_Matrix( Nmx,Nmy,Nmz , i,j+1,k );
                                AL1Dy( jL , jL+1 ) = A( Im , I_Front );
                                AL1Dy( jL , jL   ) = AL1Dy( jL , jL ) - A( Im , I_Front );
                            end
                            % Setting the value of matrix coarse nodes to 0.0
                            if ( any(CIy == j) )
                                if ~any(Well.Cells == Im)
                                    AL1Dy( jL , :   ) = 0;
                                    AL1Dy( jL , jL  ) = 1;
                                    qL1Dy( jL       ) = 0;
                                end
                            end
                            % Setting the value at the current well to 1.0 in AL1Dy
                            if any(Well.Cells == Im)
                                if w <= length(Inj)
                                    AL1Dy( jL , jL ) = AL1Dy( jL , jL ) + Well.PI * sum(Well.Mob, 2) * K(Im, 1);
                                    qL1Dy( jL , 1  ) = Well.PI * sum(Well.Mob, 2) * K(Im, 1) * (Well.p-P_init(Im));
                                else
                                    AL1Dy( jL , jL ) = AL1Dy( jL , jL ) + Well.PI * sum(Mob(Im,:), 2) * K(Im, 1);
                                    qL1Dy( jL , 1  ) = Well.PI * sum(Mob(Im,:), 2) * K(Im, 1) * (Well.p-P_init(Im));
                                end
                            end
                        end
                        % Dirichlet boundary condition at back
                        if by > 1
                            AL1Dy( 1 , : ) = 0;
                            AL1Dy( 1 , 1 ) = 1;
                        end
                        % Dirichlet boundary condition at front
                        if ey < Nmy
                            AL1Dy( NyL , :   ) = 0;
                            AL1Dy( NyL , NyL ) = 1;
                        end
                        % Solving 1D localization on AL1Dx
                        pL1Dy = AL1Dy\qL1Dy;
                        % Assigning pL1Dy for Dirichlet Boundary Condition in AL2Dz
                        for jL = 1 : NyL
                            ImL2D = Index_Matrix( NxL,NyL,1, iLcx2D(L1Dx),jL,1 );
                            AL2Dz( ImL2D , :     ) = 0;
                            AL2Dz( ImL2D , ImL2D ) = 1;
                            qL2Dz( ImL2D , 1     ) = pL1Dy(jL);
                        end
                    end
                    % Dirichlet boundary condition at left for AL2Dz
                    if bx > 1
                        for jL = 1 : NyL
                            I_Left = Index_Matrix( NxL,NyL,1 , 1,jL,1 );
                            AL2Dz( I_Left , :      ) = 0;
                            AL2Dz( I_Left , I_Left ) = 1;
                        end
                    end
                    % Dirichlet boundary condition at right for AL2Dz
                    if ex < Nmx
                        for jL = 1 : NyL
                            I_Right = Index_Matrix( NxL,NyL,1 , NxL,jL,1 );
                            AL2Dz( I_Right , :       ) = 0;
                            AL2Dz( I_Right , I_Right ) = 1;
                        end
                    end
                    % Dirichlet boundary condition at back for AL2Dz
                    if by > 1
                        for iL = 1 : NxL
                            I_Back = Index_Matrix( NxL,NyL,1 , iL,1,1 );
                            AL2Dz( I_Back , :      ) = 0;
                            AL2Dz( I_Back , I_Back ) = 1;
                        end
                    end
                    % Dirichlet boundary condition at front for AL2Dz
                    if ey < Nmy
                        for iL = 1 : NxL
                            I_Front = Index_Matrix( NxL,NyL,1 , iL,NyL,1 );
                            AL2Dz( I_Front , :       ) = 0;
                            AL2Dz( I_Front , I_Front ) = 1;
                        end
                    end
                    % Solving 2D localization on AL2Dz
                    pL2Dz = AL2Dz\qL2Dz;
                    % Assigning pL2Dz for Dirichlet Boundary Condition in AL3D
                    for jL = 1 : NyL
                        for iL = 1 : NxL
                            ImL2D  = Index_Matrix( NxL,NyL,1   , iL,jL,1            );
                            ImL3D  = Index_Matrix( NxL,NyL,NzL , iL,jL,kLcz3D(L2Dz) );
                            AL3D( ImL3D , :     ) = 0;
                            AL3D( ImL3D , ImL3D ) = 1;
                            qL3D( ImL3D , 1     ) = pL2Dz( ImL2D , 1 );
                        end
                    end 
                end
                
                %% Localization in 2D-y faces
                for L2Dy = 1:length(jLcy3D)
                    AL2Dy = zeros( NxL*NzL     );
                    qL2Dy = zeros( NxL*NzL , 1 );
                    j = by+jLcy3D(L2Dy)-1;
                    iLcx2D = [];  kLcz2D = [];
                    % AL2Dy
                    for kL = 1 : NzL
                        k = bz+kL-1;
                        for iL = 1 : NxL
                            i = bx+iL-1;
                            Im    = Index_Matrix( Nmx,Nmy,Nmz , i ,j,k  );
                            ImL2D = Index_Matrix( NxL,1  ,NzL , iL,1,kL );
                            if ( iL > 1 )
                                I_Left   = Index_Matrix( Nmx,Nmy,Nmz , i -1,j,k  );
                                I_Left_L = Index_Matrix( NxL,1  ,NzL , iL-1,1,kL );
                                AL2Dy( ImL2D , I_Left_L ) = A    ( Im    , I_Left );
                                AL2Dy( ImL2D , ImL2D    ) = AL2Dy( ImL2D , ImL2D  ) - A( Im , I_Left );
                            end
                            if ( iL < NxL )
                                I_Right   = Index_Matrix( Nmx,Nmy,Nmz , i +1,j,k  );
                                I_Right_L = Index_Matrix( NxL,1  ,NzL , iL+1,1,kL );
                                AL2Dy( ImL2D , I_Right_L ) = A    ( Im    , I_Right );
                                AL2Dy( ImL2D , ImL2D     ) = AL2Dy( ImL2D , ImL2D   ) - A( Im , I_Right );
                            end
                            if ( kL > 1 )
                                I_Bottom   = Index_Matrix( Nmx,Nmy,Nmz , i ,j,k -1 );
                                I_Bottom_L = Index_Matrix( NxL,1  ,NzL , iL,1,kL-1 );
                                AL2Dy( ImL2D , I_Bottom_L ) = A    ( Im    , I_Bottom );
                                AL2Dy( ImL2D , ImL2D      ) = AL2Dy( ImL2D , ImL2D  ) - A( Im , I_Bottom );
                            end
                            if ( kL < NzL )
                                I_Top   = Index_Matrix( Nmx,Nmy,Nmz , i ,j,k +1 );
                                I_Top_L = Index_Matrix( NxL,1  ,NzL , iL,1,kL+1 );
                                AL2Dy( ImL2D , I_Top_L ) = A    ( Im    , I_Top );
                                AL2Dy( ImL2D , ImL2D   ) = AL2Dy( ImL2D , ImL2D   ) - A( Im , I_Top );
                            end
                            % Getting the local indexes of inner coarse nodes in 2D
                            if any(CIx == i) && any(CIz == k)
                                if ( i~=bx ) && ( i~=ex ) && ( ~any(iLcx2D==iL) )
                                    iLcx2D = [iLcx2D , iL];
                                end
                                if ( k~=bz ) && ( k~=ez ) && ( ~any(kLcz2D==kL) )
                                    kLcz2D = [kLcz2D , kL];
                                end
                                % Setting the value of matrix coarse nodes to 0.0
                                if ~any(Well.Cells == Im)
                                    AL2Dy( ImL2D , :     ) = 0;
                                    AL2Dy( ImL2D , ImL2D ) = 1;
                                    AL2Dy( ImL2D         ) = 0;
                                end
                            end
                            % Setting the value at the current well to 1.0 in AL2Dy
                            if any(Well.Cells == Im)
                                if w <= length(Inj)
                                    AL2Dy( ImL2D , ImL2D ) = AL2Dy( ImL2D , ImL2D ) + Well.PI * sum(Well.Mob, 2) * K(Im, 1);
                                    qL2Dy( ImL2D , 1     ) = Well.PI * sum(Well.Mob, 2) * K(Im, 1) * (Well.p-P_init(Im));
                                else
                                    AL2Dy( ImL2D , ImL2D ) = AL2Dy( ImL2D , ImL2D ) + Well.PI * sum(Mob(Im,:), 2) * K(Im, 1);
                                    qL2Dy( ImL2D , 1     ) = Well.PI * sum(Mob(Im,:), 2) * K(Im, 1) * (Well.p-P_init(Im));
                                end
                            end
                        end
                    end
                    % Localization in 1D-x
                    for L1Dz = 1 : length(kLcz2D)
                        AL1Dx = zeros( NxL     );
                        qL1Dx = zeros( NxL , 1 );
                        k = bz+kLcz2D(L1Dz)-1;
                        % AL1Dx
                        for iL = 1 : NxL
                            i = bx+iL-1;
                            Im = Index_Matrix( Nmx,Nmy,Nmz , i,j,k );
                            if ( iL > 1 )
                                I_Left = Index_Matrix( Nmx,Nmy,Nmz , i-1,j ,k );
                                AL1Dx( iL , iL-1 ) = A( Im , I_Left );
                                AL1Dx( iL , iL   ) = AL1Dx( iL , iL ) - A( Im , I_Left );
                            end
                            if ( iL < NxL )
                                I_Right = Index_Matrix( Nmx,Nmy,Nmz , i+1,j ,k );
                                AL1Dx( iL , iL+1 ) = A( Im , I_Right );
                                AL1Dx( iL , iL   ) = AL1Dx( iL , iL ) - A( Im , I_Right );
                            end
                            % Setting the value of matrix coarse nodes to 0.0
                            if ( any(CIx == i) )
                                if ~any(Well.Cells == Im)
                                    AL1Dx( iL , :   ) = 0;
                                    AL1Dx( iL , iL  ) = 1;
                                    qL1Dx( iL       ) = 0;
                                end
                            end
                            % Setting the value at the current well to 1.0 in AL1Dx
                            if any(Well.Cells == Im)
                                if w <= length(Inj)
                                    AL1Dx( iL , iL ) = AL1Dx( iL , iL ) + Well.PI * sum(Well.Mob, 2) * K(Im, 1);
                                    qL1Dx( iL , 1  ) = Well.PI * sum(Well.Mob, 2) * K(Im, 1) * (Well.p-P_init(Im));
                                else
                                    AL1Dx( iL , iL ) = AL1Dx( iL , iL ) + Well.PI * sum(Mob(Im,:), 2) * K(Im, 1);
                                    qL1Dx( iL , 1  ) = Well.PI * sum(Mob(Im,:), 2) * K(Im, 1) * (Well.p-P_init(Im));
                                end
                            end
                        end
                        % Dirichlet boundary condition at left
                        if bx > 1
                            AL1Dx( 1 , : ) = 0;
                            AL1Dx( 1 , 1 ) = 1;
                        end
                        % Dirichlet boundary condition at right
                        if ex < Nmx
                            AL1Dx( NxL , :   ) = 0;
                            AL1Dx( NxL , NxL ) = 1;
                        end
                        % Solving 1D localization on AL1Dx
                        pL1Dx = AL1Dx\qL1Dx;
                        % Assigning pL1Dx for Dirichlet Boundary Condition in AL2Dy
                        for iL = 1 : NxL
                            ImL2D = Index_Matrix( NxL,1,NzL, iL,1,kLcz2D(L1Dz) );
                            AL2Dy( ImL2D , :     ) = 0;
                            AL2Dy( ImL2D , ImL2D ) = 1;
                            qL2Dy( ImL2D , 1     ) = pL1Dx(iL);
                        end
                    end
                    % Localization in 1D-z
                    for L1Dx = 1 : length(iLcx2D)
                        AL1Dz = zeros( NzL     );
                        qL1Dz = zeros( NzL , 1 );
                        i = bx+iLcx2D(L1Dx)-1;
                        % AL1Dz
                        for kL = 1 : NzL
                            k = bz+kL-1;
                            Im = Index_Matrix( Nmx,Nmy,Nmz , i,j,k );
                            if ( kL > 1 )
                                I_Bottom = Index_Matrix( Nmx,Nmy,Nmz , i,j,k-1 );
                                AL1Dz( kL , kL-1 ) = A( Im , I_Bottom );
                                AL1Dz( kL , kL   ) = AL1Dz( kL , kL ) - A( Im , I_Bottom );
                            end
                            if ( kL < NzL )
                                I_Top = Index_Matrix( mx,Nmy,Nmz , i,j,k+1 );
                                AL1Dz( kL , kL+1 ) = A( Im , I_Top );
                                AL1Dz( kL , kL   ) = AL1Dz( kL , kL ) - A( Im , I_Top );
                            end
                            % Setting the value of matrix coarse nodes to 0.0
                            if ( any(CIz == k) )
                                if ~any(Well.Cells == Im)
                                    AL1Dz( kL , :   ) = 0;
                                    AL1Dz( kL , kL  ) = 1;
                                    qL1Dz( kL       ) = 0;
                                end
                            end
                            % Setting the value at the current well to 1.0 in AL1Dy
                            if any(Well.Cells == Im)
                                if w <= length(Inj)
                                    AL1Dz( kL , kL ) = AL1Dz( kL , kL ) + Well.PI * sum(Well.Mob, 2) * K(Im, 1);
                                    qL1Dz( kL , 1  ) = Well.PI * sum(Well.Mob, 2) * K(Im, 1) * (Well.p-P_init(Im));
                                else
                                    AL1Dz( kL , kL ) = AL1Dz( kL , kL ) + Well.PI * sum(Mob(Im,:), 2) * K(Im, 1);
                                    qL1Dz( kL , 1  ) = Well.PI * sum(Mob(Im,:), 2) * K(Im, 1) * (Well.p-P_init(Im));
                                end
                            end
                        end
                        % Dirichlet boundary condition at bototm
                        if bz > 1
                            AL1Dz( 1 , : ) = 0;
                            AL1Dz( 1 , 1 ) = 1;
                        end
                        % Dirichlet boundary condition at top
                        if ez < Nmz
                            AL1Dz( NzL , :   ) = 0;
                            AL1Dz( NzL , NzL ) = 1;
                        end
                        % Solving 1D localization on AL1Dx
                        pL1Dz = AL1Dz\qL1Dz;
                        % Assigning pL1Dz for Dirichlet Boundary Condition in AL2Dy
                        for kL = 1 : NzL
                            ImL2D = Index_Matrix( NxL,1,NzL, iLcx2D(L1Dx),1,kL );
                            AL2Dy( ImL2D , :     ) = 0;
                            AL2Dy( ImL2D , ImL2D ) = 1;
                            qL2Dy( ImL2D , 1     ) = pL1Dz(kL);
                        end
                    end
                    % Dirichlet boundary condition at left for AL2Dy
                    if bx > 1
                        for kL = 1 : NzL
                            I_Left = Index_Matrix( NxL,1,NzL , 1,1,kL );
                            AL2Dy( I_Left , :      ) = 0;
                            AL2Dy( I_Left , I_Left ) = 1;
                        end
                    end
                    % Dirichlet boundary condition at right for AL2Dz
                    if ex < Nmx
                        for kL = 1 : NzL
                            I_Right = Index_Matrix( NxL,1,NzL , NxL,1,kL );
                            AL2Dy( I_Right , :       ) = 0;
                            AL2Dy( I_Right , I_Right ) = 1;
                        end
                    end
                    % Dirichlet boundary condition at back for AL2Dy
                    if bz > 1
                        for iL = 1 : NxL
                            I_Bottom = Index_Matrix( NxL,1,NzL , iL,1,1 );
                            AL2Dy( I_Bottom , :        ) = 0;
                            AL2Dy( I_Bottom , I_Bottom ) = 1;
                        end
                    end
                    % Dirichlet boundary condition at front for AL2Dy
                    if ez < Nmz
                        for iL = 1 : NxL
                            I_Top = Index_Matrix( NxL,1,NzL , iL,1,NzL );
                            AL2Dy( I_Top , :     ) = 0;
                            AL2Dy( I_Top , I_Top ) = 1;
                        end
                    end
                    % Solving 2D localization on AL2Dy
                    pL2Dy = AL2Dy\qL2Dy;
                    % Assigning pL2Dy for Dirichlet Boundary Condition in AL3D
                    for kL = 1 : NzL
                        for iL = 1 : NxL
                            ImL2D  = Index_Matrix( NxL,1,NzL   , iL,1,kL            );
                            ImL3D  = Index_Matrix( NxL,NyL,NzL , iL,jLcy3D(L2Dy),kL );
                            AL3D( ImL3D , :     ) = 0;
                            AL3D( ImL3D , ImL3D ) = 1;
                            qL3D( ImL3D , 1     ) = pL2Dy( ImL2D , 1 );
                        end
                    end 
                end
                
                %% Localization in 2D-x faces
                for L2Dx = 1:length(iLcx3D)
                    AL2Dx = zeros( NyL*NzL     );
                    qL2Dx = zeros( NyL*NzL , 1 );
                    i = bx+iLcx3D(L2Dx)-1;
                    jLcy2D = [];  kLcz2D = [];
                    % AL2Dx
                    for kL = 1 : NzL
                        k = bz+kL-1;
                        for jL = 1 : NyL
                            j = by+jL-1;
                            Im    = Index_Matrix( Nmx,Nmy,Nmz , i,j ,k  );
                            ImL2D = Index_Matrix( 1  ,NyL,NzL , 1,jL,kL );
                            if ( jL > 1 )
                                I_Back   = Index_Matrix( Nmx,Nmy,Nmz , i,j -1,k  );
                                I_Back_L = Index_Matrix( 1  ,NyL,NzL , 1,jL-1,kL );
                                AL2Dx( ImL2D , I_Back_L ) = A    ( Im    , I_Back );
                                AL2Dx( ImL2D , ImL2D    ) = AL2Dx( ImL2D , ImL2D  ) - A( Im , I_Back );
                            end
                            if ( jL < NyL )
                                I_Front   = Index_Matrix( Nmx,Nmy,Nmz , i,j +1,k  );
                                I_Front_L = Index_Matrix( 1  ,NyL,NzL , i,jL+1,kL );
                                AL2Dx( ImL2D , I_Front_L ) = A    ( Im    , I_Front );
                                AL2Dx( ImL2D , ImL2D     ) = AL2Dx( ImL2D , ImL2D   ) - A( Im , I_Front );
                            end
                            if ( kL > 1 )
                                I_Bottom   = Index_Matrix( Nmx,Nmy,Nmz , i,j ,k -1 );
                                I_Bottom_L = Index_Matrix( 1  ,NyL,NzL , 1,jL,kL-1 );
                                AL2Dx( ImL2D , I_Bottom_L ) = A    ( Im    , I_Bottom );
                                AL2Dx( ImL2D , ImL2D      ) = AL2Dx( ImL2D , ImL2D  ) - A( Im , I_Bottom );
                            end
                            if ( kL < NzL )
                                I_Top   = Index_Matrix( Nmx,Nmy,Nmz , i,j ,k +1 );
                                I_Top_L = Index_Matrix( 1  ,NyL,NzL , 1,jL,kL+1 );
                                AL2Dx( ImL2D , I_Top_L ) = A    ( Im    , I_Top );
                                AL2Dx( ImL2D , ImL2D   ) = AL2Dx( ImL2D , ImL2D   ) - A( Im , I_Top );
                            end
                            % Getting the local indexes of inner coarse nodes in 2D
                            if any(CIy == j) && any(CIz == k)
                                if ( j~=by ) && ( j~=ey ) && ( ~any(jLcy2D==jL) )
                                    jLcy2D = [jLcy2D , jL];
                                end
                                if ( k~=bz ) && ( k~=ez ) && ( ~any(kLcz2D==kL) )
                                    kLcz2D = [kLcz2D , kL];
                                end
                                % Setting the value of matrix coarse nodes to 0.0
                                if ~any(Well.Cells == Im)
                                    AL2Dx( ImL2D , :     ) = 0;
                                    AL2Dx( ImL2D , ImL2D ) = 1;
                                    AL2Dx( ImL2D         ) = 0;
                                end
                            end
                            % Setting the value at the current well to 1.0 in AL2Dx
                            if any(Well.Cells == Im)
                                if w <= length(Inj)
                                    AL2Dx( ImL2D , ImL2D ) = AL2Dx( ImL2D , ImL2D ) + Well.PI * sum(Well.Mob, 2) * K(Im, 1);
                                    qL2Dx( ImL2D , 1     ) = Well.PI * sum(Well.Mob, 2) * K(Im, 1) * (Well.p-P_init(Im));
                                else
                                    AL2Dx( ImL2D , ImL2D ) = AL2Dx( ImL2D , ImL2D ) + Well.PI * sum(Mob(Im,:), 2) * K(Im, 1);
                                    qL2Dx( ImL2D , 1     ) = Well.PI * sum(Mob(Im,:), 2) * K(Im, 1) * (Well.p-P_init(Im));
                                end
                            end
                        end
                    end
                    % Localization in 1D-y
                    for L1Dz = 1 : length(kLcz2D)
                        AL1Dy = zeros( NyL     );
                        qL1Dy = zeros( NyL , 1 );
                        k = bz+kLcz2D(L1Dz)-1;
                        % AL1Dy
                        for jL = 1 : NyL
                            j = by+jL-1;
                            Im = Index_Matrix( Nmx,Nmy,Nmz , i,j,k );
                            if ( jL > 1 )
                                I_back = Index_Matrix( Nmx,Nmy,Nmz , i,j-1,k );
                                AL1Dy( jL , jL-1 ) = A( Im , I_back );
                                AL1Dy( jL , jL   ) = AL1Dy( jL , jL ) - A( Im , I_back );
                            end
                            if ( jL < NyL )
                                I_Front = Index_Matrix( Nmx,Nmy,Nmz , i,j+1,k );
                                AL1Dy( jL , jL+1 ) = A( Im , I_Front );
                                AL1Dy( jL , jL   ) = AL1Dy( jL , jL ) - A( Im , I_Front );
                            end
                            % Setting the value of matrix coarse nodes to 0.0
                            if ( any(CIy == j) )
                                if ~any(Well.Cells == Im)
                                    AL1Dy( jL , :   ) = 0;
                                    AL1Dy( jL , jL  ) = 1;
                                    qL1Dy( jL       ) = 0;
                                end
                            end
                            % Setting the value at the current well to 1.0 in AL1Dy
                            if any(Well.Cells == Im)
                                if w <= length(Inj)
                                    AL1Dy( jL , jL ) = AL1Dy( jL , jL ) + Well.PI * sum(Well.Mob, 2) * K(Im, 1);
                                    qL1Dy( jL , 1  ) = Well.PI * sum(Well.Mob, 2) * K(Im, 1) * (Well.p-P_init(Im));
                                else
                                    AL1Dy( jL , jL ) = AL1Dx( jL , jL ) + Well.PI * sum(Mob(Im,:), 2) * K(Im, 1);
                                    qL1Dy( jL , 1  ) = Well.PI * sum(Mob(Im,:), 2) * K(Im, 1) * (Well.p-P_init(Im));
                                end
                            end
                        end
                        % Dirichlet boundary condition at back
                        if by > 1
                            AL1Dy( 1 , : ) = 0;
                            AL1Dy( 1 , 1 ) = 1;
                        end
                        % Dirichlet boundary condition at front
                        if ey < Nmy
                            AL1Dy( NyL , :   ) = 0;
                            AL1Dy( NyL , NyL ) = 1;
                        end
                        % Solving 1D localization on AL1Dy
                        pL1Dy = AL1Dy\qL1Dy;
                        % Assigning pL1Dx for Dirichlet Boundary Condition in AL2Dz
                        for jL = 1 : NyL
                            ImL2D = Index_Matrix( 1,NyL,NzL , 1,jL,kLcz2D(L1Dz) );
                            AL2Dx( ImL2D , :     ) = 0;
                            AL2Dx( ImL2D , ImL2D ) = 1;
                            qL2Dx( ImL2D , 1     ) = pL1Dy(jL);
                        end
                    end
                    % Localization in 1D-z
                    for L1Dy = 1 : length(jLcy2D)
                        AL1Dz = zeros( NzL     );
                        qL1Dz = zeros( NzL , 1 );
                        j = bx+jLcy2D(L1Dy)-1;
                        % AL1Dz
                        for kL = 1 : NzL
                            k = bz+kL-1;
                            Im = Index_Matrix( Nmx,Nmy,Nmz , i,j,k );
                            if ( kL > 1 )
                                I_Bottom = Index_Matrix( Nmx,Nmy,Nmz , i,j,k-1 );
                                AL1Dz( kL , kL-1 ) = A( Im , I_Bottom );
                                AL1Dz( kL , kL   ) = AL1Dz( kL , kL ) - A( Im , I_Bottom );
                            end
                            if ( kL < NzL )
                                I_Top = Index_Matrix( Nmx,Nmy,Nmz , i,j,k+1 );
                                AL1Dz( kL , kL+1 ) = A( Im , I_Top );
                                AL1Dz( kL , kL   ) = AL1Dz( kL , kL ) - A( Im , I_Top );
                            end
                            % Setting the value of matrix coarse nodes to 0.0
                            if ( any(CIz == k) )
                                if ~any(Well.Cells == Im)
                                    AL1Dz( kL , :   ) = 0;
                                    AL1Dz( kL , kL  ) = 1;
                                    qL1Dz( kL       ) = 0;
                                end
                            end
                            % Setting the value at the current well to 1.0 in AL1Dz
                            if any(Well.Cells == Im)
                                if w <= length(Inj)
                                    AL1Dz( kL , kL ) = AL1Dz( kL , kL ) + Well.PI * sum(Well.Mob, 2) * K(Im, 1);
                                    qL1Dz( kL , 1  ) = Well.PI * sum(Well.Mob, 2) * K(Im, 1) * (Well.p-P_init(Im));
                                else
                                    AL1Dz( kL , kL ) = AL1Dz( kL , kL ) + Well.PI * sum(Mob(Im,:), 2) * K(Im, 1);
                                    qL1Dz( kL , 1  ) = Well.PI * sum(Mob(Im,:), 2) * K(Im, 1) * (Well.p-P_init(Im));
                                end
                            end
                        end
                        % Dirichlet boundary condition at bottom
                        if bz > 1
                            AL1Dz( 1 , : ) = 0;
                            AL1Dz( 1 , 1 ) = 1;
                        end
                        % Dirichlet boundary condition at top
                        if ez < Nmz
                            AL1Dz( NzL , :   ) = 0;
                            AL1Dz( NzL , NzL ) = 1;
                        end
                        % Solving 1D localization on AL1Dx
                        pL1Dz = AL1Dz\qL1Dz;
                        % Assigning pL1Dz for Dirichlet Boundary Condition in AL2Dx
                        for kL = 1 : NzL
                            ImL2D = Index_Matrix( 1,NyL,NzL , 1,jLcy2D(L1Dy),kL );
                            AL2Dx( ImL2D , :     ) = 0;
                            AL2Dx( ImL2D , ImL2D ) = 1;
                            qL2Dx( ImL2D , 1     ) = pL1Dz(kL);
                        end
                    end
                    % Dirichlet boundary condition at back for AL2Dx
                    if by > 1
                        for kL = 1 : NzL
                            I_Back = Index_Matrix( 1,NyL,NzL , 1,1,kL );
                            AL2Dx( I_Back , :      ) = 0;
                            AL2Dx( I_Back , I_Back ) = 1;
                        end
                    end
                    % Dirichlet boundary condition at front for AL2Dx
                    if ey < Nmy
                        for kL = 1 : NzL
                            I_Front = Index_Matrix( 1,NyL,NzL , 1,NyL,kL );
                            AL2Dx( I_Front , :       ) = 0;
                            AL2Dx( I_Front , I_Front ) = 1;
                        end
                    end
                    % Dirichlet boundary condition at bottom for AL2Dx
                    if bz > 1
                        for jL = 1 : NyL
                            I_Bottom = Index_Matrix( 1,NyL,NzL , 1,jL,1 );
                            AL2Dx( I_Bottom , :        ) = 0;
                            AL2Dx( I_Bottom , I_Bottom ) = 1;
                        end
                    end
                    % Dirichlet boundary condition at top for AL2Dx
                    if ez < Nmz
                        for jL = 1 : NyL
                            I_Top = Index_Matrix( 1,NyL,NzL , 1,jL,NzL );
                            AL2Dx( I_Top , :     ) = 0;
                            AL2Dx( I_Top , I_Top ) = 1;
                        end
                    end
                    % Solving 2D localization on AL2Dz
                    pL2Dx = AL2Dx\qL2Dx;
                    % Assigning pL2Dx for Dirichlet Boundary Condition in AL3D
                    for kL = 1 : NzL
                        for jL = 1 : NyL
                            ImL2D = Index_Matrix( 1,NyL,NzL   , 1           ,jL,kL );
                            ImL3D = Index_Matrix( NxL,NyL,NzL , iLcx3D(L2Dx),jL,kL );
                            AL3D( ImL3D , :     ) = 0;
                            AL3D( ImL3D , ImL3D ) = 1;
                            qL3D( ImL3D , 1     ) = pL2Dx( ImL2D , 1 );
                        end
                    end 
                end
                
                %% Dirichlet boundary condition at left for AL3D
                if bx > 1
                    for kL = 1 : NzL
                        for jL = 1 : NyL
                            I_Left = Index_Matrix( NxL,NyL,NzL , 1,jL,kL );
                            AL3D( I_Left , :      ) = 0;
                            AL3D( I_Left , I_Left ) = 1;
                        end
                    end
                end
                % Dirichlet boundary condition at left for AL3D
                if ex < Nmx
                    for kL = 1 : NzL
                        for jL = 1 : NyL
                            I_Right = Index_Matrix( NxL,NyL,NzL , NxL,jL,kL );
                            AL3D( I_Right , :       ) = 0;
                            AL3D( I_Right , I_Right ) = 1;
                        end
                    end
                end
                % Dirichlet boundary condition at back for AL3D
                if by > 1
                    for kL = 1 : NzL
                        for iL = 1 : NxL
                            I_Back = Index_Matrix( NxL,NyL,NzL , iL,1,kL );
                            AL3D( I_Back , :      ) = 0;
                            AL3D( I_Back , I_Back ) = 1;
                        end
                    end
                end
                % Dirichlet boundary condition at front for AL3D
                if ey < Nmy
                    for kL = 1 : NzL
                        for iL = 1 : NxL
                            I_Front = Index_Matrix( NxL,NyL,NzL , iL,NyL,kL );
                            AL3D( I_Front , :       ) = 0;
                            AL3D( I_Front , I_Front ) = 1;
                        end
                    end
                end
                % Dirichlet boundary condition at bottom for AL3D
                if bz > 1
                    for jL = 1 : NyL
                        for iL = 1 : NxL
                            I_Bottom = Index_Matrix( NxL,NyL,NzL , iL,jL,1 );
                            AL3D( I_Bottom , :        ) = 0;
                            AL3D( I_Bottom , I_Bottom ) = 1;
                        end
                    end
                end
                % Dirichlet boundary condition at top for AL3D
                if ez < Nmz
                    for jL = 1 : NyL
                        for iL = 1 : NxL
                            I_Top = Index_Matrix( NxL,NyL,NzL , iL,jL,NzL );
                            AL3D( I_Top , :     ) = 0;
                            AL3D( I_Top , I_Top ) = 1;
                        end
                    end
                end

                %% Obtaining Solution:
                pL3D = AL3D\qL3D;
                
                % Adding solutions to pressure vector
                for kL = 1 : NzL
                    k = bz+kL-1;
                    for jL = 1 : NyL
                        j = by+jL-1;
                        for iL = 1 : NxL
                            i = bx+iL-1;
                            Im    = Index_Matrix( Nmx,Nmy,Nmz , i ,j ,k  );
                            ImL3D = Index_Matrix( NxL,NyL,NzL , iL,jL,kL );
                            deltaP_w( Im ) = pL3D( ImL3D );
                        end
                    end
                end
            end % End of loop over the wells
        end % End of StaticMultilevelPressureGuess
    end
end