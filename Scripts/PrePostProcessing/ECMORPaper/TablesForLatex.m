% Pressure and Saturation errors
DSs = {'DS001', 'DS003', 'DS005', 'DS007', 'DS01', 'DS02', 'DS03'};
Dss = [0.01; 0.03; 0.05; 0.07; 0.1; 0.2; 0.3];
PressureError = zeros(6, 8);
SaturationError = zeros(6, 8);
for dsi=1:7
    DS = DSs{dsi};
    ErrorComparison_2;
%Active nodes
    PressureError(dsi, 1) = mean(Const2L.Stats(:,4));
    PressureError(dsi, 2) = mean(Bilin2L.Stats(:,4));
    PressureError(dsi, 3) = mean(MS2L.Stats(:,4));

    SaturationError(dsi, 1) = mean(Const2L.Stats(:,4));
    SaturationError(dsi, 2) = mean(Bilin2L.Stats(:,4));
    SaturationError(dsi, 3) = mean(MS2L.Stats(:,4));
   
%Errors
    PressureError(dsi,5) = Const2L_EP;
    PressureError(dsi,6) = Bilin2L_EP;
    PressureError(dsi,7) = MS2L_EP;

    SaturationError(dsi,5) = Const2L_ES;
    SaturationError(dsi,6) = Bilin2L_ES;
    SaturationError(dsi,7) = MS2L_ES;

end
%ANvsDs
file = fopen('../../ownCloud/Matteo_PhD/Papers/ADM_JCP/ADM_review1/2016_02_01/ADM_2016_02_01/Figures/TABLES/SPE10BANvsDsC2L.dat', 'w');
fprintf(file, '%10.5f %10.5f\n', [Dss, PressureError(:,1)]');
fclose(file);
file = fopen('../../ownCloud/Matteo_PhD/Papers/ADM_JCP/ADM_review1/2016_02_01/ADM_2016_02_01/Figures/TABLES/SPE10BANvsDsB2L.dat', 'w');
fprintf(file, '%10.5f %10.5f\n', [Dss, PressureError(:,2)]');
fclose(file);
file = fopen('../../ownCloud/Matteo_PhD/Papers/ADM_JCP/ADM_review1/2016_02_01/ADM_2016_02_01/Figures/TABLES/SPE10BANvsDsMS2L.dat', 'w');
fprintf(file, '%10.2f %10.2f\n', [Dss, PressureError(:,3)]');
fclose(file);

file = fopen('../../ownCloud/Matteo_PhD/Papers/ADM_JCP/ADM_review1/2016_02_01/ADM_2016_02_01/Figures/TABLES/SPE10BEPC2L.dat', 'w');
fprintf(file, '%10.5f %10.5f\n', [Dss, PressureError(:,5)]');
fclose(file);
file = fopen('../../ownCloud/Matteo_PhD/Papers/ADM_JCP/ADM_review1/2016_02_01/ADM_2016_02_01/Figures/TABLES/SPE10BEPB2L.dat', 'w');
fprintf(file, '%10.5f %10.5f\n', [Dss, PressureError(:,6)]');
fclose(file);
file = fopen('../../ownCloud/Matteo_PhD/Papers/ADM_JCP/ADM_review1/2016_02_01/ADM_2016_02_01/Figures/TABLES/SPE10BEPMS2L.dat', 'w');
fprintf(file, '%10.5f %10.5f\n', [Dss, PressureError(:,7)]');
fclose(file);

file = fopen('../../ownCloud/Matteo_PhD/Papers/ADM_JCP/ADM_review1/2016_02_01/ADM_2016_02_01/Figures/TABLES/SPE10BESC2L.dat', 'w');
fprintf(file, '%10.5f %10.5f\n', [Dss, SaturationError(:,5)]');
fclose(file);
file = fopen('../../ownCloud/Matteo_PhD/Papers/ADM_JCP/ADM_review1/2016_02_01/ADM_2016_02_01/Figures/TABLES/SPE10BESB2L.dat', 'w');
fprintf(file, '%10.5f %10.5f\n', [Dss, SaturationError(:,6)]');
fclose(file);
file = fopen('../../ownCloud/Matteo_PhD/Papers/ADM_JCP/ADM_review1/2016_02_01/ADM_2016_02_01/Figures/TABLES/SPE10BESMS2L.dat', 'w');
fprintf(file, '%10.5f %10.5f\n', [Dss, SaturationError(:,7)]');
fclose(file);
