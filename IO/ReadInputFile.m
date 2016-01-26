%READ INPUT File
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%ReadInputFile
fileID = fopen(InputFile, 'r');
%// Read lines from input file
inputMatrix = textscan(fileID, '%s', 'Delimiter', '\n');
fclose(fileID);
temp = strfind(inputMatrix{1}, 'TITLE'); % Search a specific string and find all rows containing matches
Problem = char(inputMatrix{1}(find(~cellfun('isempty', temp)) + 1));

%%%%%%%%%%%%%PROPERTIES OF THE RESERVOIR%%%%%%%%%%%%%%%%
temp = strfind(inputMatrix{1}, 'DIMENS'); % Search a specific string and find all rows containing matches
size = find(~cellfun('isempty', temp));
temp = strfind(inputMatrix{1}, 'SPECGRID');
grid = find(~cellfun('isempty', temp));
temp = strfind(inputMatrix{1}, 'PERMX'); 
perm = find(~cellfun('isempty', temp));
temp = strfind(inputMatrix{1}, 'POR'); 
por = find(~cellfun('isempty', temp));
temp = strfind(inputMatrix{1}, 'TOTALTIME');
x = find(~cellfun('isempty', temp));
T = str2double(inputMatrix{1}(x + 1))*24*3600;
[Grid, K]=ReservoirProperties(size, grid, perm, por, inputMatrix{1});
clear temp size grid perm por x

%%%%%%%%%%%%%FLUID PROPERTIES%%%%%%%%%%%%%%%%
fplot = 0;    %if 1 the fractional function curves are plotted
temp = strfind(inputMatrix{1}, 'VISCOSITY'); % Search a specific string and find all rows containing matches
viscosity = find(~cellfun('isempty', temp));
temp = strfind(inputMatrix{1}, 'RELPERM');
relperm = find(~cellfun('isempty', temp));
temp = strfind(inputMatrix{1}, 'CAPILLARITY');
capillarity = find(~cellfun('isempty', temp));
temp = strfind(inputMatrix{1}, 'FOAM');
foam = find(~cellfun('isempty', temp));
temp = strfind(inputMatrix{1}, 'COREY');
corey = find(~cellfun('isempty', temp));
Fluid = FluidProperties(viscosity, relperm, capillarity, foam, corey, inputMatrix{1});
clear temp viscosity relperm capillarity foam corey;

%%%%%%%%%%%%%WELLS%%%%%%%%%%%%%%%%
temp = regexp(inputMatrix{1}, 'INJ\d', 'match'); 
inj = find(~cellfun('isempty', temp));
temp = regexp(inputMatrix{1}, 'PROD\d', 'match');
prod = find(~cellfun('isempty', temp));
[Inj, Prod] = WellsProperties(inj, prod, inputMatrix{1}, Grid,K);
clear temp inputMatrix;

%%%%%%%%%%%%%%%SIMULATOR'S SETTINGS%%%%%%%%%%%
InputFile = strcat(InputDirectory, '/SimulatorSettings.txt');
fileID = fopen(InputFile, 'r');
%// Read lines from input file
inputMatrix = textscan(fileID, '%s', 'Delimiter', '\n');
fclose(fileID);
temp = strfind(inputMatrix{1}, 'FIM'); % Search a specific string and find all rows containing matches
settings = find(~cellfun('isempty', temp));
if settings ~= 0
    Strategy = 'FIM';
    impsat = 0;
else
    Strategy = 'Sequential';
    temp = strfind(inputMatrix{1}, 'SEQUENTIAL'); % Search a specific string and find all rows containing matches
    settings = find(~cellfun('isempty', temp));
    temp = strfind(inputMatrix{1}, 'IMPSAT');
    impsat = find(~cellfun('isempty', temp));
end
temp = strfind(inputMatrix{1}, 'TIMESTEPS');
x = find(~cellfun('isempty', temp));
TimeSteps = str2double(inputMatrix{1}(x+1));
temp = strfind(inputMatrix{1}, 'ADM');
adm = find(~cellfun('isempty', temp));
[FIM, Sequential, ADMSettings] = ...
    SimulatorSettings(TimeSteps, Strategy, settings, impsat, adm, inputMatrix);

%%%%Properties of Injected fluid%%%%
for i=1:length(inj)
    [Inj(i).Mw, Inj(i).Mo, Inj(i).dMw, Inj(i).dMo] = Mobilities(1,Fluid);
    Inj(i).water = zeros(TimeSteps+1,1);
end
for i=1:length(prod)
    Prod(i).water = zeros(TimeSteps+1,1);
    Prod(i).oil = zeros(TimeSteps+1,1);
end

%%%%%%%%%%%%%OPTIONS%%%%%%%%%%%%%%%%
temp = strfind(inputMatrix{1}, 'OUTPUT');
x = find(~cellfun('isempty', temp));
Options.PlotSolution = char(inputMatrix{1}(x+1)); %Matlab or VTK
Options.Pressure_3D = 0; % 0 or 1, if 1 pressure plot in 3D
Options.problem_1D = 0; % if 1, the plotting for a 1D problem is activated
Options.ContourPlot = 0; % If 1 dynamic contour plot
Options.PlotResiduals = 0; % If 1 Residuals are plotted

clear settings impsat adm inputMatrix x