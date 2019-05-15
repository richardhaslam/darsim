%% pEDFM-ADM 2D 135x135 Homogeneous with 30 MixedConductiveFractures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % ADM tol 05
% Directory = 'D:\SURFdrive\Simulation\DARSim2\InputDesktopWork\pEDFM_ADM\SinglePhaseGeo_2D_135x135_30frac_pEDFM\MixedConductiveFractures\Homogeneous\ADM_05';
% Result = DARSim2ResSim(Directory,'SinglePhaseGeo.txt','D:\SURFdrive\Simulation\DARSim2\Permeability');
% Result.Simulation.DiscretizationModel.OperatorsHandler = [];
% save(strcat(Directory,'\Output\Result.mat'),'Result');
% clear all;

%% pEDFM-ADM 2D 136x136 Homogeneous with 30 MixedConductiveFractures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADM tol 05
Directory = 'D:\SURFdrive\Simulation\DARSim2\InputDesktopWork\pEDFM_ADM\SinglePhaseGeo_2D_136x136_30frac_pEDFM\MixedConductiveFractures\Homogeneous\ADM_05';
Result = DARSim2ResSim(Directory,'SinglePhaseGeo.txt','D:\SURFdrive\Simulation\DARSim2\Permeability');
Result.Simulation.DiscretizationModel.OperatorsHandler = [];
save(strcat(Directory,'\Output\Result.mat'),'Result');
clear all;
% 
% % ADM tol 10
% Directory = 'D:\SURFdrive\Simulation\DARSim2\InputDesktopWork\pEDFM_ADM\SinglePhaseGeo_2D_136x136_30frac_pEDFM\MixedConductiveFractures\Homogeneous\ADM_10';
% Result = DARSim2ResSim(Directory,'SinglePhaseGeo.txt','D:\SURFdrive\Simulation\DARSim2\Permeability');
% Result.Simulation.DiscretizationModel.OperatorsHandler = [];
% save(strcat(Directory,'\Output\Result.mat'),'Result');
% clear all;
% 
% % ADM tol 20
% Directory = 'D:\SURFdrive\Simulation\DARSim2\InputDesktopWork\pEDFM_ADM\SinglePhaseGeo_2D_136x136_30frac_pEDFM\MixedConductiveFractures\Homogeneous\ADM_20';
% Result = DARSim2ResSim(Directory,'SinglePhaseGeo.txt','D:\SURFdrive\Simulation\DARSim2\Permeability');
% Result.Simulation.DiscretizationModel.OperatorsHandler = [];
% save(strcat(Directory,'\Output\Result.mat'),'Result');
% clear all;
% 
% % ADM tol 50
% Directory = 'D:\SURFdrive\Simulation\DARSim2\InputDesktopWork\pEDFM_ADM\SinglePhaseGeo_2D_136x136_30frac_pEDFM\MixedConductiveFractures\Homogeneous\ADM_50';
% Result = DARSim2ResSim(Directory,'SinglePhaseGeo.txt','D:\SURFdrive\Simulation\DARSim2\Permeability');
% Result.Simulation.DiscretizationModel.OperatorsHandler = [];
% save(strcat(Directory,'\Output\Result.mat'),'Result');
% clear all;
% 
% % FineScale
% Directory = 'D:\SURFdrive\Simulation\DARSim2\InputDesktopWork\pEDFM_ADM\SinglePhaseGeo_2D_136x136_30frac_pEDFM\MixedConductiveFractures\Homogeneous\FineScale';
% Result = DARSim2ResSim(Directory,'SinglePhaseGeo.txt','D:\SURFdrive\Simulation\DARSim2\Permeability');
% save(strcat(Directory,'\Output\Result.mat'),'Result');
% clear all;


%% pEDFM-ADM 2D 136x136 Heterogeneous with 30 MixedConductiveFractures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % ADM tol 50
% Directory = 'D:\SURFdrive\Simulation\DARSim2\InputDesktopWork\pEDFM_ADM\SinglePhaseGeo_2D_136x136_30frac_pEDFM\MixedConductiveFractures\Heterogeneous\ADM_50';
% Result = DARSim2ResSim(Directory,'SinglePhaseGeo.txt','D:\SURFdrive\Simulation\DARSim2\Permeability');
% Result.Simulation.DiscretizationModel.OperatorsHandler = [];
% save(strcat(Directory,'\Output\Result.mat'),'Result');
% clear all;
% 
% % ADM tol 20
% Directory = 'D:\SURFdrive\Simulation\DARSim2\InputDesktopWork\pEDFM_ADM\SinglePhaseGeo_2D_136x136_30frac_pEDFM\MixedConductiveFractures\Heterogeneous\ADM_20';
% Result = DARSim2ResSim(Directory,'SinglePhaseGeo.txt','D:\SURFdrive\Simulation\DARSim2\Permeability');
% Result.Simulation.DiscretizationModel.OperatorsHandler = [];
% save(strcat(Directory,'\Output\Result.mat'),'Result');
% clear all;
% 
% % ADM tol 10
% Directory = 'D:\SURFdrive\Simulation\DARSim2\InputDesktopWork\pEDFM_ADM\SinglePhaseGeo_2D_136x136_30frac_pEDFM\MixedConductiveFractures\Heterogeneous\ADM_10';
% Result = DARSim2ResSim(Directory,'SinglePhaseGeo.txt','D:\SURFdrive\Simulation\DARSim2\Permeability');
% Result.Simulation.DiscretizationModel.OperatorsHandler = [];
% save(strcat(Directory,'\Output\Result.mat'),'Result');
% clear all;
% 
% % ADM tol 05
% Directory = 'D:\SURFdrive\Simulation\DARSim2\InputDesktopWork\pEDFM_ADM\SinglePhaseGeo_2D_136x136_30frac_pEDFM\MixedConductiveFractures\Heterogeneous\ADM_05';
% Result = DARSim2ResSim(Directory,'SinglePhaseGeo.txt','D:\SURFdrive\Simulation\DARSim2\Permeability');
% Result.Simulation.DiscretizationModel.OperatorsHandler = [];
% save(strcat(Directory,'\Output\Result.mat'),'Result');
% clear all;

% % FineScale
% Directory = 'D:\SURFdrive\Simulation\DARSim2\InputDesktopWork\pEDFM_ADM\SinglePhaseGeo_2D_136x136_30frac_pEDFM\MixedConductiveFractures\Heterogeneous\FineScale';
% Result = DARSim2ResSim(Directory,'SinglePhaseGeo.txt','D:\SURFdrive\Simulation\DARSim2\Permeability');
% save(strcat(Directory,'\Output\Result.mat'),'Result');
% clear all;