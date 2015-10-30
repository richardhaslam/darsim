%MsFV Restriction and Prolongation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MsR, MsP, C] = MSFVOperators(FineGrid, CoarseGrid, Ap, level)
%MSFV Restriction and Prolongation Operators
Nf = FineGrid.Nx*FineGrid.Ny;
Nc = CoarseGrid.Nx*CoarseGrid.Ny;

%Permutation Matrix
[G, Ni, Ne, Nn] = PermutationMatrix(FineGrid, Nf, CoarseGrid);

%Restriction and Prolongation for Fine Grid ordering
MsR = MsRestriction(FineGrid, CoarseGrid, Nf, Nc, level);
[MsP, C] = MsProlongation(Ni, Ne, Nn, Ap, G, FineGrid);
MsP = G'*MsP;
C = G'*C*G;
end

%% Permutation Matrix
function [P, nii, nee, nnn] = PermutationMatrix(FineGrid, Nf, CoarseGrid)
nxf  = FineGrid.Nx;
nyf  = FineGrid.Ny;
nf   = Nf;
nxcf = CoarseGrid.CoarseFactor(1)/FineGrid.CoarseFactor(1);
nycf = CoarseGrid.CoarseFactor(2)/FineGrid.CoarseFactor(2);
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

% if nxc == 1, 
%     nxd  = 1; 
%     jcen = 0;
% else 
%     nxd  = nxc+1;
%     jcen = ceil(nxcf/2);
% end
% if nyc == 1, 
%     nyd  = 1;
%     icen = 0;
% else 
%     nyd  = nyc+1; 
%     icen = ceil(nycf/2);
% end

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

%% MSFV Restriction Operator
function MsR = MsRestriction(FineGrid, CoarseGrid, Nf, Nc, level)
MsR = zeros(Nc, Nf);
if (FineGrid.CoarseFactor(1) == 1)
    for r=1:Nc
        %coordinates of fine cells contained in the coarse block
        Imin = CoarseGrid.I(r) - floor((CoarseGrid.CoarseFactor(1) - 1)/2);
        Imax = CoarseGrid.I(r) + ceil((CoarseGrid.CoarseFactor(1) - 1)/2);
        Jmin = CoarseGrid.J(r) - floor((CoarseGrid.CoarseFactor(2) - 1)/2);
        Jmax = CoarseGrid.J(r) + ceil((CoarseGrid.CoarseFactor(2) - 1)/2);
        i=Imin:Imax;
        j=Jmin:Jmax;
        [p,q] = meshgrid(i, j);
        pairs = [p(:), q(:)];
        %indexes of the fine cells
        c = pairs(:,1) + (pairs(:,2)-1)*FineGrid.Nx;
        %I make 1 those columns
        MsR(r,c) = 1;
    end
else
    for r=1:Nc
        for c=1:Nf
            if (r == FineGrid.Father(c, level))
                %I make 1 those columns
                MsR(r,c) = 1;
            end
        end
    end
end
MsR = sparse(MsR);
end

%% MS Prolongation Operator
function [MsP, C] = MsProlongation(Ni, Ne, Nn, Ap, G, FineGrid)
% 1. Reorder finescale system
%A_tpfa = MakeTPFA(Ap, FineGrid);
tildeA   = G*Ap*G';
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

% 9. Define matrix: when applied to rhs gives correction functions            
C = [Mii_inv        -Mii_inv*Mie*Mee_inv      sparse(Ni,Nn);        ...
     sparse(Ne,Ni)               Mee_inv      sparse(Ne,Nn);        ...
     sparse(Nn,Ni)         sparse(Nn,Ne)      sparse(Nn,Nn)];       
end

% function B = MakeTPFA(A, FineGrid)
% Nx = FineGrid.Nx;
% N = FineGrid.Nx*FineGrid.Ny;
% B=zeros(N);
% for i=1:N
%     for j=1:N
%         if (j==i || j==i-1 || j==i+1 || j==i-Nx || j==i+Nx)
%             B(i,j) = A(i,j);
%         end
%     end  
% end
% B = sparse(B);
% end
