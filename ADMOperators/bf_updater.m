%%  Basis functions updater
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef bf_updater < handle
    properties
        A
    end
    methods 
        function [MsP, MsC] = MsProlongation(obj, FineGrid, CoarseGrid, Dimensions)
            cf = CoarseGrid.CoarseFactor ./ FineGrid.CoarseFactor;
            % Permutation Matrix
            [G, Ni, Nf, Ne, Nv] = obj.PermutationMatrix(FineGrid, CoarseGrid, cf);
            % Reorder A based on dual coarse grid partition 
            tildeA = G * obj.A * G';
            [MsP, MsC] = obj.ComputeMsP(tildeA, Ni, Nf, Ne, Nv, Dimensions);
            MsP = G' * MsP;
            if ~isempty(MsC)
                MsC = G' * MsC * G;
            end
        end
        function [MsP, MsC] = ComputeMsP(obj, tildeA, Ni, Nf, Ne, Nv, Dimensions)
            switch(Dimensions)
                case (1)
                    % 2. Define edge-edge (ee) block
                    Mee = tildeA(Ni+Nf+1:Ni+Nf+Ne,Ni+Nf+1:Ni+Nf+Ne);
                    % 3. Define edge-vertex (ev) block
                    Mev = tildeA(Ni+Nf+1:Ni+Nf+Ne,Ni+Nf+Ne+1:Ni+Nf+Ne+Nv);
                    % 4. Compute inverse of (ii), (ff) and (ee) blocks
                    Mvv = speye(Nv,Nv);
                    % 1D
                    Edges = -slvblk(Mee, Mev);
                    % Prolongation Operator
                    MsP = [Edges;...
                           Mvv];
                    % Correction functions operator
                    MsC = [];
                    
                case (2)
                    % 3. Define face-x (ff,fe,fv) blocks
                    Mff = tildeA(Ni+1:Ni+Nf, Ni+1:Ni+Nf); 
                    Mfe = tildeA(Ni+1:Ni+Nf, Ni+Nf+1:Ni+Nf+Ne);
                    Mfv = tildeA(Ni+1:Ni+Nf, Ni+Nf+Ne+1:Ni+Nf+Ne+Nv);
                    % 4. Define edge-x (ee,ev) blocks
                    Mee = tildeA(Ni+Nf+1:Ni+Nf+Ne,Ni+Nf+1:Ni+Nf+Ne) + diag(sum(Mfe,1));
                    Mev = tildeA(Ni+Nf+1:Ni+Nf+Ne,Ni+Nf+Ne+1:Ni+Nf+Ne+Nv);
                    % 5. Define vertex-vertex (vv) block
                    Mvv = speye(Nv,Nv);
                    % 6. Compute inverse of (ii), (ff) and (ee) blocks
                    Edges = -slvblk(Mee, Mev);
                    Faces = -slvblk(Mff, Mfe * Edges + Mfv);
                    % 7. Obtaining Prolongation
                    MsP = [Faces;...
                           Edges;...
                           Mvv];
                    % 8. Correction functions operator
                    MsC = [];
                    
                case(3)
                    % 2. Define interior-x (ii,if,ie,iv) blocks
                    Mii = tildeA(1:Ni, 1:Ni);
                    Mif = tildeA(1:Ni, Ni+1:Ni+Nf);
                    Mie = tildeA(1:Ni, Ni+Nf+1:Ni+Nf+Ne);
                    Miv = tildeA(1:Ni, Ni+Nf+Ne+1:Ni+Nf+Ne+Nv);
                    % 3. Define face-x (ff,fe,fv) blocks
                    Mff = tildeA(Ni+1:Ni+Nf, Ni+1:Ni+Nf) + diag(sum(Mif,1));
                    Mfe = tildeA(Ni+1:Ni+Nf, Ni+Nf+1:Ni+Nf+Ne);
                    Mfv = tildeA(Ni+1:Ni+Nf, Ni+Nf+Ne+1:Ni+Nf+Ne+Nv);
                    % 4. Define edge-x (ee,ev) blocks
                    Mee = tildeA(Ni+Nf+1:Ni+Nf+Ne,Ni+Nf+1:Ni+Nf+Ne) + diag(sum(Mie,1)) + diag(sum(Mfe,1));
                    Mev = tildeA(Ni+Nf+1:Ni+Nf+Ne,Ni+Nf+Ne+1:Ni+Nf+Ne+Nv);
                    % 5. Define vertex-vertex (vv) block
                    Mvv = speye(Nv,Nv);
                    % 6. Compute inverse of (ii), (ff) and (ee) blocks
                    Edges = -slvblk(Mee, Mev);
                    Faces = -slvblk(Mff, Mfe * Edges + Mfv);
                    Interiors = -slvblk(Mii, Mif * Faces + Mie * Edges + Miv);
                    % 7. Obtaining Prolongation
                    MsP = [Interiors;...
                           Faces    ;...
                           Edges    ;...
                           Mvv     ];
                    % 8. Correction functions
                    MsC = [];
            end
        end
        function [P, nii, nff, nee, nvv] = PermutationMatrix(obj, FineGrid, CoarseGrid, cf)
            % The ordering is done based on the dual coarse scale internal data
            %
            %     Non overlappig dual grid
            %
            %                 v eeee              Dual grid 6 x 6
            %                 e iiii             we consider the first left 5 x 5 in numbering
            %                 e iiii
            %                 e iiii
            %
            %     1  2  [1]   5  6  .....          7  8  [3]   13  14
            %     3  4  [2]   9 10  ....          11 12  [4]   15  16
            %    [5][6] (1)  [7] [8]           [9] [10]  (2)   [15] [16]
            %    17 18  [11]  25 26             27 28   [17]   41   42
            %    19 20  [12]  29  30            31 32   [18]   43    44
            %
            %     iyc = 1 ixc =1                       iyc = 1   ixc =2
            %
            %    21  22  [13]  33 34  .....     35 36   [19]  45 46
            %    23  24  [14]  37 38  ...       39 40   [20]  47  4850
            % [21][22]   (3) [ 23][24]        [25] [26]  (4)   [xx7] [xx8]
            %            [27]                            [yy3]
            %            [28]                            [yy4]
            %
            %     iyc = 2 ixc =1                       iyc = 2   ixc =2
            
            %% 0. Define local variables
            nxf = vertcat(FineGrid.Nx);
            nyf = vertcat(FineGrid.Ny);
            nzf = vertcat(FineGrid.Nz);
            nf  = vertcat(FineGrid.N );
            if FineGrid(1).CoarseLevel > 0
                nxf( [FineGrid.hasCoarseNodes] == 0 ) = 0;
                nyf( [FineGrid.hasCoarseNodes] == 0 ) = 0;
                nzf( [FineGrid.hasCoarseNodes] == 0 ) = 0;
                nf ( [FineGrid.hasCoarseNodes] == 0 ) = 0;
            end
            
            nxcf = cf(:,1);
            nycf = cf(:,2);
            nzcf = cf(:,3);
            
            nxc = vertcat(CoarseGrid.Nx);  nxc( [CoarseGrid.hasCoarseNodes] == 0 ) = 0;
            nyc = vertcat(CoarseGrid.Ny);  nyc( [CoarseGrid.hasCoarseNodes] == 0 ) = 0;
            nzc = vertcat(CoarseGrid.Nz);  nzc( [CoarseGrid.hasCoarseNodes] == 0 ) = 0;
            
            %% 1. Define number of cells of each block
            if ~CoarseGrid(1).Vertex_On_Corner
                nxc_hybrid = nxc;
                nyc_hybrid = nyc;
                nzc_hybrid = nzc;
            else
                nxc_hybrid = max((nxc-1),1);
                nyc_hybrid = max((nyc-1),1);
                nzc_hybrid = max((nzc-1),1);
            end
            
            nii = 0; % interiors
            nff = 0; % faces
            nee = 0; % edges
            nvv = 0; % verteces
            
            % interiors (only in matrix, no interiors in fractures)
            nii = ((nxcf(1)-1) * (nycf(1)-1) * (nzcf(1)-1)) * nxc_hybrid(1)*nyc_hybrid(1)*nzc_hybrid(1);

            for m=1:length(FineGrid)
                % faces
                if CoarseGrid(m).hasCoarseNodes == 0 && FineGrid(m).Ny > 1
                    nff = nff + nf(m);
                else
                    nff = nff + (nxcf(m)-1)*(nycf(m)-1)*nxc_hybrid(m)*nyc_hybrid(m)*nzc(m) + ...
                              (nycf(m)-1)*(nzcf(m)-1)*nyc_hybrid(m)*nzc_hybrid(m)*nxc(m) + ...
                              (nxcf(m)-1)*(nzcf(m)-1)*nxc_hybrid(m)*nzc_hybrid(m)*nyc(m);
                end
                        
                % edges
                if CoarseGrid(m).hasCoarseNodes == 0 && FineGrid(m).Ny == 1
                    nee = nee + nf(m);
                else
                    nee = nee + (nxcf(m)-1)*nxc_hybrid(m)*nyc(m)*nzc(m) + ...
                                (nycf(m)-1)*nyc_hybrid(m)*nxc(m)*nzc(m) + ...
                                (nzcf(m)-1)*nzc_hybrid(m)*nxc(m)*nyc(m);
                end
                        
                % verteces
                nvv = nvv + nxc(m)*nyc(m)*nzc(m);
            end
            P = sparse(sum(nf),sum(nf));
            
            %% 2. Construct matrix and vector
            % Coarsening factors are assumed to be odd numbers. 
            iii0 = 0;
            iff0 = nii;
            iee0 = nii + nff;
            ivv0 = nii + nff + nee;
            
            if ~CoarseGrid(1).Vertex_On_Corner
                %% If the coarse grids are constructed normally (no coarse
                % nodes on the corners)
                nxd = nxc + 1;
                nyd = nyc + 1;
                nzd = nzc + 1;
                icen = ceil(nxcf/2);
                jcen = ceil(nycf/2);
                kcen = ceil(nzcf/2);
                % Now fill in entries of permutation matrix
                for m = 1 : length(FineGrid)
                    for k = 1:nzd(m) % k index of dual domain
                        for j = 1:nyd(m) % j index of dual domain
                            for i = 1:nxd(m) % i index of dual domain
                                for kz= 1:nzcf(m) % k index within a primal coarse cell
                                    for jy= 1:nycf(m) % j index within a primal coarse cell
                                        for ix = 1:nxcf(m) % i index within a primal coarse cell
                                            % calculate the cell index in the original matrix
                                            ii = (i-1) * nxcf(m) - icen(m) + ix;
                                            jj = (j-1) * nycf(m) - jcen(m) + jy;
                                            kk = (k-1) * nzcf(m) - kcen(m) + kz;
                                            if ( ii >= 1 && ii <= nxf(m) && jj >= 1 && jj <= nyf(m) && kk >=1 && kk<=nzf(m))
                                                % There are dual domains that are
                                                % not within the domain
                                                ijk = sum(nf(1:m-1)) + ii + (jj-1)*nxf(m) + (kk-1)*nxf(m)*nyf(m); % global index
                                                if (ix  ~= 1 && jy ~= 1 && kz ~=1)
                                                    % internal point
                                                    iii0 = iii0 + 1;
                                                    P(iii0, ijk) = 1;
                                                elseif  (ix == 1 && jy == 1 && kz==1)
                                                    % node points
                                                    ivv0 = ivv0 + 1;
                                                    P(ivv0, ijk) = 1;
                                                elseif (ix == 1 && jy == 1 && (kz ~=1 || kz~=nzcf(m)))
                                                    % edge xy points
                                                    iee0 = iee0 + 1;
                                                    P(iee0, ijk) = 1;
                                                elseif (ix == 1 && kz == 1 && (jy ~=1 || jy~=nycf(m)))
                                                    % edge xz points
                                                    iee0 = iee0 + 1;
                                                    P(iee0, ijk) = 1;
                                                elseif (jy == 1 && kz == 1 && (ix ~=1 || ix~=nxcf(m)))
                                                    % edge yz points
                                                    iee0 = iee0 + 1;
                                                    P(iee0, ijk) = 1;
                                                else
                                                    % face points
                                                    iff0 = iff0 + 1;
                                                    P(iff0, ijk) = 1;
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                    % In case there are no coarse node in this media, all
                    % the cells are faces (if 2D) or edges (if 1D)
                    if CoarseGrid(m).hasCoarseNodes == 0 && FineGrid(m).Ny > 1
                        for j = 1 : nyf(m)
                            for i = 1 : nxf(m)
                                ijk = sum(nf(1:m-1)) + (j-1) * FineGrid(m).Ny + i;
                                iff0 = iff0 + 1;
                                P(iff0, ijk) = 1;
                            end
                        end
                    end
                    if CoarseGrid(m).hasCoarseNodes == 0 && FineGrid(m).Ny == 1
                        for i = 1 : nxf(m)
                            ijk = sum(nf(1:m-1)) + i;
                            iee0 = iee0 + 1;
                            P(iee0, ijk) = 1;
                        end
                    end
                end
            else
                %% If the coarse nodes are on the corners
                icen = ceil(nxcf/2);
                jcen = ceil(nycf/2);
                kcen = ceil(nzcf/2);
                % Now fill in entries of permutation matrix
                for m = 1 : length(FineGrid)
                    for k = 1:nzc(m)
                        for j = 1:nyc(m)
                            for i = 1:nxc(m)
                                for kz= 1:nzcf(m)
                                    for jy= 1:nycf(m)
                                        for ix = 1:nxcf(m)
                                            % calculate the cell index in the original matrix
                                            ii = (i-1) * nxcf(m) - icen(m) + ix + 1;
                                            jj = (j-1) * nycf(m) - jcen(m) + jy + 1;
                                            kk = (k-1) * nzcf(m) - kcen(m) + kz + 1;
                                            if ( ii >= 1 && ii <= nxf(m) && jj >= 1 && jj <= nyf(m) && kk >=1 && kk<=nzf(m))
                                                % There are dual domains that are
                                                % not within the domain
                                                ijk = sum(nf(1:m-1)) + ii + (jj-1)*nxf(m) + (kk-1)*nxf(m)*nyf(m); % global index
                                                if (ix  ~= icen(m) && jy ~= jcen(m) && kz ~=kcen(m))
                                                    % internal point
                                                    iii0 = iii0 + 1;
                                                    P(iii0, ijk) = 1;
                                                elseif (ix == icen(m) && jy == jcen(m) && kz == kcen(m))
                                                    % node points
                                                    ivv0 = ivv0 + 1;
                                                    P(ivv0, ijk) = 1;
                                                elseif (ix == icen(m) && jy == jcen(m) && kz ~=kcen(m))
                                                    % edge xy points
                                                    iee0 = iee0 + 1;
                                                    P(iee0, ijk) = 1;
                                                elseif (ix == icen(m) && kz == kcen(m) && jy ~= jcen(m))
                                                    % edge xz points
                                                    iee0 = iee0 + 1;
                                                    P(iee0, ijk) = 1;
                                                elseif (jy == jcen(m) && kz == kcen(m) && ix ~= icen(m))
                                                    % edge yz points
                                                    iee0 = iee0 + 1;
                                                    P(iee0, ijk) = 1;
                                                else
                                                    % face points
                                                    iff0 = iff0 + 1;
                                                    P(iff0, ijk) = 1;
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                    % In case there are no coarse node in this media, all
                    % the cells are faces (if 2D) or edges (if 1D)
                    if CoarseGrid(m).hasCoarseNodes == 0 && FineGrid(m).Ny > 1
                        for j = 1 : nyf(m)
                            for i = 1 : nxf(m)
                                ijk = sum(nf(1:m-1)) + (j-1) * nyf(m) + i;
                                iff0 = iff0 + 1;
                                P(iff0, ijk) = 1;
                            end
                        end
                    end
                    if CoarseGrid(m).hasCoarseNodes == 0 && FineGrid(m).Ny == 1
                        for i = 1 : nxf(m)
                            ijk = sum(nf(1:m-1)) + i;
                            iee0 = iee0 + 1;
                            P(iee0, ijk) = 1;
                        end
                    end
                end
            end
        end
        function UpdatePressureMatrix(obj, P, Grid)
            obj.A = P' * obj.A * P;
            obj.A = obj.TransformIntoTPFA(obj.A, Grid);
        end
        function Ac = TransformIntoTPFA(obj, Ac, Grid)
            Nx = Grid.Nx;
            Ny = Grid.Ny;
            N = Grid.N;
            x1 = [0; diag(Ac, 1)];
            x2 = [diag(Ac, -1); 0];
            diagy1 = diag(Ac, Nx);
            diagy2 = diag(Ac, -Nx);
            y1 = zeros(N, 1);
            y2 = zeros(N, 1);
            y1(N - length(diagy1) + 1:end) = diagy1;
            y2(1:length(diagy1)) = diagy2;
            diagz1 = diag(Ac, Nx*Ny);
            diagz2 = diag(Ac, -Nx*Ny);
            z1 = zeros(N, 1);
            z2 = zeros(N, 1);
            z1(N - length(diagz1) + 1:end) = diagz1;
            z2(1:length(diagz1)) = diagz2;
            DiagVecs = [z2, y2, x2, -z2-y2-x2-z1-y1-x1, x1, y1, z1];
            DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
            Ac = spdiags(DiagVecs,DiagIndx,N,N);
        end
    end
end