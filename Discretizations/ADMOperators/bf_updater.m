%  Basis functions updater
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 16 August 2016
%Last modified: 16 August 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef bf_updater < handle
    properties
        A
    end
    methods 
        function MsP = MsProlongation(obj, FineGrid, CoarseGrid, cf)
            %Permutation Matrix
            [G, Ni, Nf, Ne, Nn] = obj.PermutationMatrix(FineGrid, CoarseGrid, cf);
            % 1. Reorder finescale system
            tildeA   = G * obj.A * G';
            % 2. Define interior-interior (ii) block
            Mii = tildeA(1:Ni, 1:Ni);
            % 3. define interior-face (if) block
            Mif = tildeA(1:Ni, Ni+1:Ni+Nf);
            % 4. Define face-face block
            Mff = tildeA(Ni+1:Ni+Nf, Ni+1:Ni+Nf) + diag(sum(Mif,1));
            % Define edge-face (fe) block
            Mfe = tildeA(Ni+1:Ni+Nf, Ni+Nf+1:Ni+Nf+Ne);
            % 5. Define edge-edge (ee) block
            Mee = tildeA(Ni+Nf+1:Ni+Nf+Ne,Ni+Nf+1:Ni+Nf+Ne) + diag(sum(Mfe,1));
            % 6. Define edge-node (en) block
            Men = tildeA(Ni+Nf+1:Ni+Nf+Ne,Ni+Nf+Ne+1:Ni+Nf+Ne+Nn);
            % 7. Compute inverse of (ii), (ff) and (ee) blocks
            Mii_inv = Mii^-1;      
            Mff_inv = Mff^-1;
            Mee_inv = Mee^-1;
            
            % 8. Assemble Prolongation operator: columns are the basis functions
            switch (Ni)
                case(0)
                    % 2D
                    % when 2D there are no interiors. 
                    MsP = [Mff_inv*(Mfe*Mee_inv*Men);...
                          -Mee_inv*Men;...
                          speye(Nn,Nn)];
                otherwise
                    % 3D
                    MsP = [-Mii_inv*(Mif*Mff_inv*Mfe*Mee_inv*Men);...
                          Mff_inv*(Mfe*Mee_inv*Men);...
                          -Mee_inv*Men;...
                          speye(Nn,Nn)];
            end
            MsP = G'*MsP;
        end
        %% Permutation Matrix
        function [P, nii, nff, nee, nnn] = PermutationMatrix(obj, FineGrid, CoarseGrid, cf)
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
            nxf = FineGrid.Nx;
            nyf = FineGrid.Ny;
            nzf = FineGrid.Nz;
            nf  = FineGrid.N;
            nxcf = cf(1);
            nycf = cf(2);
            nzcf = cf(3);
            
            nxc = CoarseGrid.Nx;
            nyc = CoarseGrid.Ny;
            nzc = CoarseGrid.Nz;
            
            %% 1. Define number of cells of each block
            nii =  ((nxcf-1) * (nycf-1) * (nzcf-1)) * (nxc * nyc * nzc);% interiors
            nff = ((nxcf-1)*(nycf-1) + (nycf-1)*(nzcf-1) + (nxcf-1)*(nzcf-1)) * nxc * nyc * nzc; % faces
            nee = ((nxcf-1) + (nycf-1) + (nzcf-1)) * nxc * nyc * nzc; % edges
            nnn = nxc *nyc*nzc; % verteces
            
            P = sparse(nf);
            
            %% 2. Construct matrix and vector
            % Coarsening factors are assumed to be odd numbers. 
            iii0 = 0;
            iff0 = nii;
            iee0 = nii + nff;
            inn0 = nii + nff + nee;
            
            nxd = nxc + 1;
            nyd = nyc + 1;
            nzd = nzc + 1;
            icen = ceil(nxcf/2);
            jcen = ceil(nycf/2);
            kcen = ceil(nzcf/2);
            % Now fill in entries of permutation matrix 
            for k = 1:nzd
                for j = 1:nyd
                    for i = 1:nxd
                        for kz= 1:nzcf
                            for jy= 1:nycf
                                for ix = 1:nxcf
                                    % calculate the cell index in the original matrix
                                    ii = (i-1) * nxcf - icen + ix;
                                    jj = (j-1) * nycf - jcen + jy;
                                    kk = (k-1) * nzcf - kcen + kz;
                                    if ( ii >= 1 && ii <= nxf && jj >= 1 && jj <= nyf && kk >=1 && kk<=nzf)
                                        % There are dual domains that are
                                        % not withing the domain
                                        ijk = ii + (jj-1)*nxf + (kk-1)*nxf*nyf; % global index
                                        if (ix  ~= 1 && jy ~= 1 && kz ~=1)
                                            % internal point
                                            iii0 = iii0 + 1;
                                            P(iii0, ijk) = 1;  
                                        elseif  (ix == 1 && jy == 1 && kz==1)
                                            % node points
                                            inn0 = inn0 + 1;
                                            P(inn0, ijk) = 1;
                                        elseif (ix == 1 && jy == 1 && (kz ~=1 || kz~=nzcf))
                                            % edge xy points
                                            iee0 = iee0 + 1;
                                            P(iee0, ijk) = 1;
                                        elseif (ix == 1 && kz == 1 && (jy ~=1 || jy~=nycf))
                                            % edge xz points
                                            iee0 = iee0 + 1;
                                            P(iee0, ijk) = 1;
                                        elseif (jy == 1 && kz == 1 && (ix ~=1 || ix~=nxcf))
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
        end
        function TransformIntoTPFA(obj, Nx, Ny)
            [N, ~] = size(obj.A);
            x1 = [0; diag(obj.A, 1)];
            x2 = [diag(obj.A, -1); 0];
            diagy1 = diag(obj.A, Nx);
            diagy2 = diag(obj.A, -Nx);
            y1 = zeros(N, 1);
            y2 = zeros(N, 1);
            y1(N - length(diagy1) + 1:end) = diagy1;
            y2(1:length(diagy1)) = diagy2;
            diagz1 = diag(obj.A, Nx*Ny);
            diagz2 = diag(obj.A, -Nx*Ny);
            z1 = zeros(N, 1);
            z2 = zeros(N, 1);
            z1(N - length(diagz1) + 1:end) = diagz1;
            z2(1:length(diagz1)) = diagz2;
            DiagVecs = [z2, y2,x2, -z2-y2-x2--z1-y1-x1,x1,y1,z1];
            DiagIndx = [-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny];
            obj.A = spdiags(DiagVecs,DiagIndx,N,N);
        end
    end
end
