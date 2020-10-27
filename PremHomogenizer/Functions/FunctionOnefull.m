% Permeability Homogenizer Function for DARSim Reservoir Simulator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Manuela Bastidas
% Created: 
% Last modified:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of the effective permeability associated to the problem
%                    - div(K\grad p) = f
%
% Inputs: Field - Fine-scale permeability field (correspond tyo the
%                 complete fine scale domain)
%         gridX and gridY - Coarse scale cartesiad grid points in X and Y.
%                           This allows the construction of non-uniforms
%                           coarse grids. Include 0 and L
%                           Interpretation: Position in the fine-scale
%                           matrix.
% Outputs: Effective_fine - Return the values of the effective permeability
%                           with the same resolution of the initial
%                           fine-scale grid.
%          Effective_ind - Coarse scale effective permeability. One
%                          effective tensor for each coarse-scale grid point.
%
% Remarks: This function use non-dimensional quantities
% 
% Auxiliar functions:  edge, PreProcess_Identification,
% Position_indicators, MicroSolver_FEM, EfectivePermTensor
%-------------------------------------------------------
function PermEffective = FunctionOnefull(K_perm,gridX,gridY)

longX = diff(gridX)';
longY = diff(gridY)';

PermEffective = zeros(length(gridY)-1,length(gridX)-1,1);
PermEffectiveFull = zeros(size(K_perm,1),size(K_perm,2),1);

for ii = 1:length(gridX)-1
    for jj = 1:length(gridY)-1
        %% MICRO GRID
        % Inside onf the loop to allow non-uniform meshes 
        micro_gridX = sum(longX(1:ii))-longX(ii)+1:sum(longX(1:ii));
        micro_gridY = sum(longY(1:jj))-longY(jj)+1:sum(longY(1:jj));
        [xx_micro,yy_micro] = meshgrid(micro_gridX,micro_gridY);
        
        % Create micro-grid
        Micro_geo            = struct();
        Micro_geo.coordinate = [xx_micro(:),yy_micro(:)];
        Micro_geo.element    = delaunay(xx_micro,yy_micro);
        % nElement -> number of element at each mesh
        Micro_geo.nElement   = size(Micro_geo.element,1);
        % nnodes -> number of nodes at each mesh
        Micro_geo.nnodes     = size(Micro_geo.coordinate,1);
        
        % Pre-process the micro-mesh
        [~,Micro_geo.nodes2edge,Micro_geo.noedges,Micro_geo.edge2element,Micro_geo.interioredge] = ...
            edge(Micro_geo.element,Micro_geo.coordinate);
        Micro_geo = PreProcess_Identification(Micro_geo);
        % Micro_geo = Position_indicators(Micro_geo,0);
        
        %% Fine-scale Permeability        
        K_micro = K_perm(micro_gridY,micro_gridX);
        
        %% Micro cell problems
        % NUMERICAL GRADIENT
        [Kx,Ky] = gradient(K_micro);
        
        [Micro_geo,Vel1,Vel2]    = MicroSolver_FEM(Micro_geo,xx_micro,yy_micro,K_micro,Kx,Ky);
        
        [K_eff] = EfectivePermTensor(Micro_geo,xx_micro,yy_micro,K_micro,Vel1,Vel2);
%         PermEffective(jj,ii,:) = [K_eff(1,1) K_eff(2,2)];
        PermEffective(jj,ii,:) = abs(K_eff(1,1));
        
%         PermEffectiveFull(micro_gridY,micro_gridX,1) = kron(PermEffective(jj,ii,1),ones(longX(ii),longY(jj)));
%         PermEffectiveFull(micro_gridY,micro_gridX,2) = kron(PermEffective(jj,ii,2),ones(longX(ii),longY(jj)));
    end
end