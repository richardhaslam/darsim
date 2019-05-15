% Folders Molar ADM
clc
% ADM_dfdx = '../Input/SATINTERPOLATOR/Homogeneous_1D/Small_dt/Homo1D_adm_dfdx/Output/Solution/Homo1D_dfdx_smalldt_Sol';
% ADM_dfdt = '../Input/SATINTERPOLATOR/Homogeneous_1D/Small_dt/Homo1D_adm_dfdt/Output/Solution/Homo1D_dfdt_smalldt_Sol';
% FS = '../Input/SATINTERPOLATOR/Homogeneous_1D/Small_dt/Homo1D_fs/Output/Solution/Homo1D_fs_smalldt_Sol';
% MaxTimeStep = 203;

% ADM_dfdx = '../Input/SATINTERPOLATOR/Homogeneous_2D/Small_dt/Homo2D_adm_dfdx/Output/Solution/Homo2D_dfdx_smalldt_Sol';
% ADM_dfdt = '../Input/SATINTERPOLATOR/Homogeneous_2D/Small_dt/Homo2D_adm_dfdt/Output/Solution/Homo2D_dfdt_smalldt_Sol';
% FS = '../Input/SATINTERPOLATOR/Homogeneous_2D/Small_dt/Homo2D_fs/Output/Solution/Homo2D_fs_smalldt_Sol';
% MaxTimeStep = 207;

ADM_dfdx = '../Input/SATINTERPOLATOR/SPE10B/SPE10B_dfdx/Output/Solution/SPE10B_dfdx_Sol';
ADM_dfdt = '../Input/SATINTERPOLATOR/SPE10B/SPE10B_dfdt/Output/Solution/SPE10B_dfdt_Sol';
FS = '../Input/SATINTERPOLATOR/SPE10B/SPE10B_fs/Output/Solution/SPE10B_fs_Sol';
MaxTimeStep = 40;

E_P_dfdx = zeros(MaxTimeStep, 1);
E_S_dfdx = zeros(MaxTimeStep, 1);

E_P_dfdt = zeros(MaxTimeStep, 1);
E_S_dfdt = zeros(MaxTimeStep, 1);

for time = 1:MaxTimeStep
    % Error norm for time t
    FS_sol = load(strcat(FS, num2str(time),'.txt'));
    ADM_dfdx_sol = load(strcat(ADM_dfdx, num2str(time),'.txt'));
    ADM_dfdt_sol = load(strcat(ADM_dfdt, num2str(time),'.txt'));
    
    E_P_dfdx(time) = norm(FS_sol(:,2) - ADM_dfdx_sol(:,2), 2)/norm(FS_sol(:,2), 2);
    E_S_dfdx(time) = norm(FS_sol(:,3) - ADM_dfdx_sol(:,3), 2)/norm(FS_sol(:,3), 2);
    E_P_dfdt(time) = norm(FS_sol(:,2) - ADM_dfdt_sol(:,2), 2)/norm(FS_sol(:,2), 2);
    E_S_dfdt(time) = norm(FS_sol(:,3) - ADM_dfdt_sol(:,3), 2)/norm(FS_sol(:,3), 2);
end
% Average over times
P_error_dfdx = mean(E_P_dfdx)
S_error_dfdx = mean(E_S_dfdx)

P_error_dfdt = mean(E_P_dfdt)
S_error_dfdt = mean(E_S_dfdt)


