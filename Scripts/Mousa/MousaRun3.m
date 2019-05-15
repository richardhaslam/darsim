%% pEDFM-ADM 2D 136x136 Homogeneous with 30 MixedConductiveFractures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % FS
% Directory = 'D:\SURFdrive\Simulation\DARSim2\InputDesktopWork\pEDFM_ADM\ImmHomo_2D_136x136_30frac_pEDFM\MixedConductiveFractures\Homogeneous\FineScale';
% Result = DARSim2ResSim(Directory,'ImmHomo.txt','D:\SURFdrive\Simulation\DARSim2\Permeability');
% save(strcat(Directory,'\Output\Result.mat'),'Result');
% clear all;

% ADM tol 0.1
Directory = 'D:\SURFdrive\Simulation\DARSim2\InputDesktopWork\pEDFM_ADM\ImmHomo_2D_136x136_30frac_pEDFM\MixedConductiveFractures\Homogeneous\ADM_0.1';
Result = DARSim2ResSim(Directory,'ImmHomo.txt','D:\SURFdrive\Simulation\DARSim2\Permeability');
save(strcat(Directory,'\Output\Result.mat'),'Result');
clear all;

% ADM tol 0.3
Directory = 'D:\SURFdrive\Simulation\DARSim2\InputDesktopWork\pEDFM_ADM\ImmHomo_2D_136x136_30frac_pEDFM\MixedConductiveFractures\Homogeneous\ADM_0.3';
Result = DARSim2ResSim(Directory,'ImmHomo.txt','D:\SURFdrive\Simulation\DARSim2\Permeability');
save(strcat(Directory,'\Output\Result.mat'),'Result');
clear all;

% ADM tol 0.5
Directory = 'D:\SURFdrive\Simulation\DARSim2\InputDesktopWork\pEDFM_ADM\ImmHomo_2D_136x136_30frac_pEDFM\MixedConductiveFractures\Homogeneous\ADM_0.5';
Result = DARSim2ResSim(Directory,'ImmHomo.txt','D:\SURFdrive\Simulation\DARSim2\Permeability');
save(strcat(Directory,'\Output\Result.mat'),'Result');
clear all;

% ADM tol 0.8
Directory = 'D:\SURFdrive\Simulation\DARSim2\InputDesktopWork\pEDFM_ADM\ImmHomo_2D_136x136_30frac_pEDFM\MixedConductiveFractures\Homogeneous\ADM_0.8';
Result = DARSim2ResSim(Directory,'ImmHomo.txt','D:\SURFdrive\Simulation\DARSim2\Permeability');
save(strcat(Directory,'\Output\Result.mat'),'Result');
clear all;


%% pEDFM-ADM 2D 136x136 Heterogeneous with 30 MixedConductiveFractures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FS
Directory = 'D:\SURFdrive\Simulation\DARSim2\InputDesktopWork\pEDFM_ADM\ImmHomo_2D_136x136_30frac_pEDFM\MixedConductiveFractures\Heterogeneous\FineScale';
Result = DARSim2ResSim(Directory,'ImmHomo.txt','D:\SURFdrive\Simulation\DARSim2\Permeability');
save(strcat(Directory,'\Output\Result.mat'),'Result');
clear all;

% ADM tol 0.1
Directory = 'D:\SURFdrive\Simulation\DARSim2\InputDesktopWork\pEDFM_ADM\ImmHomo_2D_136x136_30frac_pEDFM\MixedConductiveFractures\Heterogeneous\ADM_0.1';
Result = DARSim2ResSim(Directory,'ImmHomo.txt','D:\SURFdrive\Simulation\DARSim2\Permeability');
save(strcat(Directory,'\Output\Result.mat'),'Result');
clear all;

% ADM tol 0.3
Directory = 'D:\SURFdrive\Simulation\DARSim2\InputDesktopWork\pEDFM_ADM\ImmHomo_2D_136x136_30frac_pEDFM\MixedConductiveFractures\Heterogeneous\ADM_0.3';
Result = DARSim2ResSim(Directory,'ImmHomo.txt','D:\SURFdrive\Simulation\DARSim2\Permeability');
save(strcat(Directory,'\Output\Result.mat'),'Result');
clear all;

% ADM tol 0.5
Directory = 'D:\SURFdrive\Simulation\DARSim2\InputDesktopWork\pEDFM_ADM\ImmHomo_2D_136x136_30frac_pEDFM\MixedConductiveFractures\Heterogeneous\ADM_0.5';
Result = DARSim2ResSim(Directory,'ImmHomo.txt','D:\SURFdrive\Simulation\DARSim2\Permeability');
save(strcat(Directory,'\Output\Result.mat'),'Result');
clear all;

% ADM tol 0.8
Directory = 'D:\SURFdrive\Simulation\DARSim2\InputDesktopWork\pEDFM_ADM\ImmHomo_2D_136x136_30frac_pEDFM\MixedConductiveFractures\Heterogeneous\ADM_0.8';
Result = DARSim2ResSim(Directory,'ImmHomo.txt','D:\SURFdrive\Simulation\DARSim2\Permeability');
save(strcat(Directory,'\Output\Result.mat'),'Result');
clear all;