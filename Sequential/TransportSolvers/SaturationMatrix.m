%Upwind Flux matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = SaturationMatrix(Grid, U, q)
%Builds Upwind Flux matrix
Nx = Grid.Nx; 
Ny = Grid.Ny; 
N = Grid.N;                                    % number of unknowns
fp = min(q,0);                                % production >> it's a negative flux
% − flow in negative coordinate - direction (XN,YN) 
XN = min(U.x,0); x1 = reshape(XN(1:Nx,:),N,1);  
YN = min(U.y,0); y1 = reshape(YN(:,1:Ny),N,1);

% − flow in positive coordinate (XP, YP)
XP = max(U.x,0); x2 = reshape(XP(2:Nx+1,:),N,1);
YP = max(U.y,0); y2 = reshape(YP(:,2:Ny+1),N,1);

DiagVecs = [y2,x2,fp+x1-x2+y1-y2,-x1,-y1]; % diagonal vectors
DiagIndx = [-Nx,-1,0,1,Nx]; % diagonal index
A = spdiags(DiagVecs,DiagIndx,N,N);      
end