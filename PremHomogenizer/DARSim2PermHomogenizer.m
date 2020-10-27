% Permeability Homogenizer Function for DARSim Reservoir Simulator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Manuela Bastidas & Mousa HosseiniMehr
% Created: 
% Last modified: 27 October 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Perm_CoarseScale_Homogenized = DARSim2PermHomogenizer(N_FineScale, MaxCoarseningLevel, CoarseningRatio, Perm_FineScale)
%% Run homogenization
TotalStart = tic;

% Adding the homogenized permeability for each coarsening level only for matrix
disp(newline)
disp('------------------------')
disp('Homogenization of permeability for DLGR:');

% By now it only works in 2D
Perm_CoarseScale_Homogenized = cell(MaxCoarseningLevel,1);
for CL = 1:MaxCoarseningLevel
    % HOMOGENIZATION
    fprintf('Calculating effective permeability - coarsening level %i ... ',CL);
%     N_CoarseScale = N_FineScale ./ (CoarseningRatio.^CL);
%     if any(mod(N_CoarseScale,1))
%         N_CoarseScale = (N_FineScale - 1) ./ (CoarseningRatio.^CL) + 1;
%     end
    
    gridX = linspace( 0 , N_FineScale(1) , N_FineScale(1)/(CoarseningRatio(1)^CL)+1 );
    gridY = linspace( 0 , N_FineScale(2) , N_FineScale(2)/(CoarseningRatio(2)^CL)+1 );
    K_temp = reshape( Perm_FineScale(:,1) , N_FineScale(1) , N_FineScale(2) ) .* 1e15;
    K_temp = FunctionOnefull( K_temp' , gridX , gridY );
    K_temp = K_temp';
    Perm_CoarseScale_Homogenized{CL} = K_temp(:) .* 1e-15;
    fprintf('Done\n');
end

TotalTime = toc(TotalStart);
disp(['The total time spent for permeability homogenization is ' num2str(TotalTime) ' [s].']);
disp('------------------------')
end