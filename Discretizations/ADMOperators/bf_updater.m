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
            [G, Ni, Ne, Nn] = obj.PermutationMatrix(FineGrid, CoarseGrid, cf);

            % 1. Reorder finescale system
            tildeA   = G * obj.A * G';
            %tildeA_tpfa = G*A_tpfa*G';
            % 2. Define interior-interior (ii) block
            Mii = tildeA(1:Ni, 1:Ni);
            % 3. define interior-edge (ie) block
            Mie = tildeA(1:Ni, Ni+1:Ni+Ne);
            % 4.Define interior-node (in) block
            Min = tildeA(1:Ni, Ni+Ne+1:end);
            % 5. Define edge-edge (ee) block
            Mee = tildeA(Ni+1:Ni+Ne,Ni+1:Ni+Ne) + diag(sum(Mie,1));
            %Mee = RemoveOffDiag(Mee);
            % 6. Define edge-node (en) block
            Men = tildeA(Ni+1:Ni+Ne,Ni+Ne+1:Ni+Ne+Nn);
            % 7. Compute inverse of (ii) and (ee) blocks
            Mii_inv = Mii^-1;
            Mee_inv = Mee^-1;
            
            % 8. Assemble Prolongation operator: columns are the basis functions
            MsP = [Mii_inv*(Mie*Mee_inv*Men-Min);                                       ...
                -Mee_inv*Men;                                       ...
                speye(Nn,Nn)];
        end
        %% Permutation Matrix
        function [P, nii, nee, nnn] = PermutationMatrix(obj, FineGrid, CoarseGrid, cf)
            nxf  = FineGrid.Nx;
            nyf  = FineGrid.Ny;
            nf   = FineGrid.N;
            nxcf = cf(1);
            nycf = cf(2);
            nxc  = CoarseGrid.Nx;
            nyc  = CoarseGrid.Ny;
            
            nii = (nxcf-1)*(nycf-1)*nxc*nyc;
            nee = ((nxcf-1) + (nycf -1)) * nxc *nyc;
            nnn = nxc *nyc;
            
            % DUALPERMUTATIONOPERATOR Permutation operator for reordering based on the dual grid
            %
            %  [P,nii,nee,nnn] = dualPermutationOperator(Nc,NfperC,type)
            %
            % -------------------------------------------------------------------------
            % TYPE = 1 (Default option)
            % The ordering is done based on the dual coarse scale internal data
            %
            %     Non overlappig dual grid
            %
            %                 v eeee              Dual grid 6 x 6
            %                 e iiii             we consider the first left 5 x 5 in numbering
            %                 e iiii
            %                 e iiii
            %
            %
            %     1  2  [1]   5  6  .....          7  8  [3]   13  14
            %     3  4  [2]   9 10  ....          11 12  [4]   15  16
            %    [5][6] (1)  [7] [8]           [9] [10]  (2)   [15] [16]
            %    17 18  [11]  25 26             27 28   [17]   41   42
            %    19 20  [12]  29  30            31 32   [18]   43    44
            %
            %     iyc = 1 ixc =1                       iyc = 1   ixc =2
            %
            %
            %    21  22  [13]  33 34  .....     35 36   [19]  45 46
            %    23  24  [14]  37 38  ...       39 40   [20]  47  4850
            % [21][22]   (3) [ 23][24]        [25] [26]  (4)   [xx7] [xx8]
            %            [27]                            [yy3]
            %            [28]                            [yy4]
            %
            %     iyc = 2 ixc =1                       iyc = 2   ixc =2
            %
            
            P = sparse(nf,nf);
            
            %% Construct matrix and vector
            % !!! No flow boundary condition everywhere !!!
            %   we also assume nxcf and nycf are odd numbers
            %
            
            iii0 = 0;
            iee0 = nii;
            inn0 = nii+nee;

            nyd  = nyc+1;
            nxd  = nxc+1;
            icen = ceil(nycf/2);
            jcen = ceil(nxcf/2);
            for i = 1:nyd
                for j = 1:nxd
                    for iy= 1:nycf
                        for jx= 1:nxcf
                            
                            %  calculate  the cell index in the original matrix
                            ii=  (i-1) *nycf  -icen + iy;
                            jj=  (j-1)* nxcf  -jcen + jx;
                            if ( ii >=1 && ii <=nyf && jj >=1 && jj <=nxf)
                                ij = (ii-1)*nxf + jj;
                                %%%             internal point
                                if ( iy  ~= 1 && jx ~= 1)
                                    iii0= iii0+1;
                                    P(iii0,ij ) =1;
                                    %%%            node points
                                elseif  ( iy == 1 && jx == 1)
                                    inn0=inn0+1;
                                    P(inn0, ij) =1;
                                    %%%            edge points
                                else
                                    iee0=iee0+1;
                                    P(iee0, ij) =1;
                                end
                            end
                        end
                    end
                end
            end
        end

    end
end