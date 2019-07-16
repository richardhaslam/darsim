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
disp(['Reading input file ', InputFile, ' from ', Directory]);
disp(char(5));

[Grid_N,maxLevel,Coarsening,K_original] = reader_Homogenizer(Directory, InputFile, PermDir);

%% Run homogenization
TotalStart = tic;

% Adding the homogenized permeability for each coarsening level
% only for matrix
disp(newline)
disp('------------------------')
disp('Homogenization: Calculating effective permeability');
% By now it only works in 2D
for c = 1:maxLevel
    % HOMOGENIZATION
    gridX = linspace(0,Grid_N(1),Grid_N(1)/(Coarsening(1)^c)+1);
    gridY = linspace(0,Grid_N(2),Grid_N(2)/(Coarsening(2)^c)+1);
    K_temp = reshape(K_original(:,1),Grid_N(1),Grid_N(2));
    [K_temp2] = FunctionOnefull(K_temp',gridX,gridY);
    K_temp2 = K_temp2';
    Homogenized_Perm{c} = K_temp2(:);
    fprintf('Effective permeability - coarsening level %i\n',c)
end

TotalTime = toc(TotalStart);

%% Display elapsed time
disp('------------------------')
disp(['The total time for Homogenization is ' num2str(TotalTime) ' s']); fprintf('\n');

%% Output Results

for c = 1:maxLevel
    name =sprintf('Perm_Coarse_L%i.txt',c);
    File = strcat(PermDir,'/',name);
    delete(File);
    fid = fopen(File,'a+');
    fprintf(fid, '%1.6e\n',Grid_N(1)/Coarsening(1)^c);
    fprintf(fid, '%1.6e\n',Grid_N(2)/Coarsening(2)^c);
    fprintf(fid, '%1.6e\n',Grid_N(3)/Coarsening(3)^c);
    fprintf(fid, '%1.6e\n',Homogenized_Perm{c});
    fclose('all');
end
end