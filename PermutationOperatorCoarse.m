function Gc = PermutationOperatorCoarse(CoarseGrid, FineGrid)
% Permutation operator based on coarse grid ordering
 
% The ordering is done based on the primal coarse grid  
%                  
%
%     1  2  3  4  5    .....          26 27 28 29 30
%     6  7  8  9 10    ....           31 32 33 34 35 
%    11 12[13]14 15                   36 37[38]39 40 
%    16 17 18 19 20                   41 42 43 44 45
%    21 22 23 24 25                   46 47 48 49 50
%
%     iyc = 1 ixc =1               iyc = 1   ixc =2   
   
nxf  = FineGrid.Nx;
nyf  = FineGrid.Ny;
nf   = nxf * nyf;
nxcf = CoarseGrid.CoarseFactor(1);
nycf = CoarseGrid.CoarseFactor(2);
nxc  = CoarseGrid.Nx;
nyc  = CoarseGrid.Ny;
nxycf= nxcf *nycf; 
Gc = zeros(nf,nf);

%% Construct matrix 
% !!! No flow boundary condition everywhere !!!

for i = 1:nyc
    for j = 1:nxc
        for iy= 1:nycf
            for jx= 1:nxcf         
%  calculate  the cell index in the original matrix
                ijc= ((i-1)* nxc +j-1)*nxycf + (iy-1)*nxcf+ jx;  
                ii=  (i-1)* nycf + iy ; 
                jj=  (j-1)* nxcf  + jx;
                ijf = (ii-1)*nxf + jj;
                Gc(ijc, ijf) = 1;
            end
        end
    end
end
Gc = sparse(Gc);
end