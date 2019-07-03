% Permeability Homogenizer Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: 
%Created: 
%Last modified: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DARSim2PermHomogenizer(Directory, InputFile, PermDir)
close all; clc;

% Remove some warnings 
warning('off', 'MATLAB:singularMatrix');
warning('off', 'MATLAB:nearlySingularMatrix');

%% Print title
disp('********************************************************************************');
disp('*************** Permeability Homogenizer for DARSim2 Simulator *****************');
disp('********************************************************************************');
disp(char(5));
disp(['Reading input file ', File, ' from ', Directory]);
disp(char(5));

%% Read Input File

% Build objects

% Print info to screen

%% Run pEDFM Generator
TotalStart = tic;

% Adding the homogenized permeability for each coarsening level
% only for matrix
Km = ProductionSystem.Reservoir.K*1e+15;
disp(newline)
disp('------------------------')
disp('Homogenization: Calculating effective permeability');
for c = 1:obj.maxLevel(1)
    % HOMOGENIZATION
    gridX = linspace(0,obj.Coarsening(1,1,c)*obj.CoarseGrid(1,c).Nx,obj.CoarseGrid(1,c).Nx+1);
    gridY = linspace(0,obj.Coarsening(1,2,c)*obj.CoarseGrid(1,c).Ny,obj.CoarseGrid(1,c).Ny+1);
    K_temp = reshape(Km(:,1),obj.Coarsening(1,1,c)*obj.CoarseGrid(1,c).Nx,obj.Coarsening(1,2,c)*obj.CoarseGrid(1,c).Ny);
    [K_temp2,~] = FunctionOnefull(K_temp',gridX,gridY);
    K_temp2 = K_temp2';
    obj.CoarseGrid(1,c).Perm = repmat(K_temp2(:),1,3)*1e-15;
    fprintf('Effective permeability - coarsening level %i\n',c)
end

TotalTime = toc(TotalStart);

%% Display elapsed time
disp('------------------------')
disp(['The total computation time is ' num2str(TotalTime) ' s']); fprintf('\n');

%% Output Results


end