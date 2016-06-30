%READ INPUT File
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [Problem, T, Grid, K, Fluid, Inj, Prod, Strategy, FIM, Sequential, ADMSettings, FlashSettings, TimeSteps, Options]...
    = ReadInputFile(InputDirectory, InputFile)
global Errors

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
temp = strfind(inputMatrix{1}, 'PERTURB'); 
pert = find(~cellfun('isempty', temp));
temp = strfind(inputMatrix{1}, 'POR'); 
por = find(~cellfun('isempty', temp));
temp = strfind(inputMatrix{1}, 'TEMPERATURE (K)'); 
temperature = find(~cellfun('isempty', temp));
temp = strfind(inputMatrix{1}, 'TOTALTIME');
xv = find(~cellfun('isempty', temp));
T = str2double(inputMatrix{1}(xv + 1))*24*3600;
[Grid, K] = ReservoirProperties(size, grid, perm, pert, por, temperature, inputMatrix{1});
clear temp size grid perm por temperature x

%%%%%%%%%%%%%FLUID PROPERTIES%%%%%%%%%%%%%%%%
temp = strfind(inputMatrix{1}, 'DENSITY (kg/m^3)'); % Search a specific string and find all rows containing matches
density = find(~cellfun('isempty', temp));
temp = strfind(inputMatrix{1}, 'VISCOSITY (Pa sec)'); 
viscosity = find(~cellfun('isempty', temp));
temp = strfind(inputMatrix{1}, 'RELPERM');
relperm = find(~cellfun('isempty', temp));
temp = strfind(inputMatrix{1}, 'COMPRESSIBILITY (1/Pa)');
compressibility = find(~cellfun('isempty', temp));
temp = strfind(inputMatrix{1}, 'CAPILLARITY');
capillarity = find(~cellfun('isempty', temp));
temp = strfind(inputMatrix{1}, 'FOAM');
foam = find(~cellfun('isempty', temp));
temp = strfind(inputMatrix{1}, 'COMPOSITION TYPE');
Comp_Type = find(~cellfun('isempty', temp));
temp = strfind(inputMatrix{1}, 'COMPONENT PROPERTIES');
Comp_Prop = find(~cellfun('isempty', temp));
Fluid = FluidProperties(density, viscosity, relperm, compressibility, capillarity, foam, Comp_Type, Comp_Prop, inputMatrix{1});
clear temp density viscosity relperm compressibility capillarity foam corey Comp_Type Comp_Prop;

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
xv = find(~cellfun('isempty', temp));
TimeSteps = str2double(inputMatrix{1}(xv+1));
temp = strfind(inputMatrix{1}, 'ADM');
adm = find(~cellfun('isempty', temp));
[FIM, Sequential, ADMSettings,FlashSettings] = ...
    SimulatorSettings(TimeSteps, Strategy, settings, impsat, adm, inputMatrix);

CheckInputErrors;

%%%%Properties of Injected fluid%%%%
if Errors == 0
    
    for i=1:length(inj)
        Inj(i).z = [1 0];
        switch (Inj(i).type)
            case('PressureConstrained')
                [Inj(i).rho, ~] = LinearDensity(Inj(i).p, Fluid.c, Fluid.rho);
                if (strcmp(Fluid.Type,'Immiscible') == 1)
                    Inj(i).x1 = [1 0];
                else
                    [Inj(i), SinglePhase] = Flash(Inj(i), Fluid, Grid.Tres, FlashSettings.TolFlash);
                    [Inj(i)] = ComputePhaseSaturation(Inj(i), SinglePhase);
                end
            case('RateConstrained')
                if (strcmp(Fluid.Type,'Immiscible')==1)
                    Inj(i).x1 = [1 0];
                else
                    display('Sorry, no rate constrained injectors for compositional problems.')
                end
        end
        Inj(i).x2 = 1 - Inj(i).x1;
        [Inj(i).Mw, Inj(i).Mo, Inj(i).dMw, Inj(i).dMo] = Mobilities(Inj(i).s ,Fluid);
    end
    
end

%%%%%%%%%%%%%OPTIONS%%%%%%%%%%%%%%%%
temp = strfind(inputMatrix{1}, 'OUTPUT');
xv = find(~cellfun('isempty', temp));
Options.PlotSolution = char(inputMatrix{1}(xv+1)); %Matlab or VTK


% Plot denisty range
Pdummy = linspace(min(Prod(1).p), max(Inj(1).p), 100);
[Rho, ~] = LinearDensity( Pdummy, Fluid.c, Fluid.rho);
figure(10)
plot(Pdummy, Rho(:,1), Pdummy, Rho(:,2))
xlabel('Pressure [Pa]')
ylabel('Density [kg/m^3]')
title('Density Range')
legend('Phase 1', 'Phase 2')
axis([min(Prod(1).p), max(Inj(1).p), min(min(Rho))*.9, max(max(Rho))*1.1])

clear settings impsat adm inputMatrix xv
end