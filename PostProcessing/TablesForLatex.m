% Pressure and Saturation errors
DSs = {'DS001', 'DS003', 'DS005', 'DS007', 'DS01', 'DS02', 'DS03'};
Dss = [0.01; 0.03; 0.05; 0.07; 0.1; 0.2; 0.3];
PressureError = zeros(7, 8);
SaturationError = zeros(7, 8);
for dsi=1:7
    DS = DSs{dsi};
    ErrorComparison_2;
    %Active nodes
    PressureError(dsi, 1) = mean(sum(Const2L.Stats(:,4:6), 2)./FineCells*100);
    PressureError(dsi, 2) = mean(sum(Bilin2L.Stats(:,4:6), 2)./FineCells*100);
    PressureError(dsi, 3) = mean(sum(MS2L.Stats(:,4:6), 2)./FineCells*100);
%     PressureError(dsi, 4) = mean(Bilin3L.Stats(:,4));
    SaturationError(dsi, 1) = mean(Const2L.Stats(:,4));
    SaturationError(dsi, 2) = mean(Bilin2L.Stats(:,4));
    SaturationError(dsi, 3) = mean(MS2L.Stats(:,4));
%     SaturationError(dsi, 4) = mean(Bilin3L.Stats(:,4));
    %Errors
    PressureError(dsi,5) = Const2L_EP;
    PressureError(dsi,6) = Bilin2L_EP;
    PressureError(dsi,7) = MS2L_EP;
%     PressureError(dsi,8) = Bilin3L_EP;
    SaturationError(dsi,5) = Const2L_ES;
    SaturationError(dsi,6) = Bilin2L_ES;
    SaturationError(dsi,7) = MS2L_ES;
%     SaturationError(dsi,8) = Bilin3L_ES;
end
%ANvsDs
file = fopen('../../Papers/ADM/Figures/TABLES/SPE10TANvsDsC2L.dat', 'w');
fprintf(file, '%10.5f %10.5f\n', [Dss, PressureError(:,1)]');
fclose(file);
file = fopen('../../Papers/ADM/Figures/TABLES/SPE10TANvsDsB2L.dat', 'w');
fprintf(file, '%10.5f %10.5f\n', [Dss, PressureError(:,2)]');
fclose(file);
file = fopen('../../Papers/ADM/Figures/TABLES/SPE10TANvsDsMS2L.dat', 'w');
fprintf(file, '%10.2f %10.2f\n', [Dss, PressureError(:,3)]');
fclose(file);
% file = fopen('../../Papers/ADM/Figures/TABLES/SPE10TANvsDsB3L.dat', 'w');
% fprintf(file, '%10.2f %10.2f\n', [Dss, PressureError(:,4)]');
% fclose(file);

file = fopen('../../Papers/ADM/Figures/TABLES/SPE10TEPC2L.dat', 'w');
fprintf(file, '%10.5f %10.5f\n', [Dss, PressureError(:,5)]');
fclose(file);
file = fopen('../../Papers/ADM/Figures/TABLES/SPE10TEPB2L.dat', 'w');
fprintf(file, '%10.5f %10.5f\n', [Dss, PressureError(:,6)]');
fclose(file);
file = fopen('../../Papers/ADM/Figures/TABLES/SPE10TEPMS2L.dat', 'w');
fprintf(file, '%10.5f %10.5f\n', [Dss, PressureError(:,7)]');
fclose(file);
% file = fopen('../../Papers/ADM/Figures/TABLES/SPE10TEPB3L.dat', 'w');
% fprintf(file, '%10.2f %10.2f\n', [PressureError(:,4), PressureError(:,8)]');
% fclose(file);
file = fopen('../../Papers/ADM/Figures/TABLES/SPE10TESC2L.dat', 'w');
fprintf(file, '%10.5f %10.5f\n', [Dss, SaturationError(:,5)]');
fclose(file);
file = fopen('../../Papers/ADM/Figures/TABLES/SPE10TESB2L.dat', 'w');
fprintf(file, '%10.5f %10.5f\n', [Dss, SaturationError(:,6)]');
fclose(file);
file = fopen('../../Papers/ADM/Figures/TABLES/SPE10TESMS2L.dat', 'w');
fprintf(file, '%10.5f %10.5f\n', [Dss, SaturationError(:,7)]');
fclose(file);
% file = fopen('../../Papers/ADM/Figures/TABLES/SPE10TESB3L.dat', 'w');
% fprintf(file, '%10.2f %10.2f\n', [SaturationError(:,4), SaturationError(:,8)]');
% fclose(file);