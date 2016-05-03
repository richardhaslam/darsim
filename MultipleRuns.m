%%%%Multiple Runs%%%%
disp('Start multiple runs');
disp('Fine scale');
InputDirectory = '../Input/GeoStat';
InputFile = strcat(InputDirectory, '/GeoStat.txt');
ResSimulator;

clear all
disp('ADM simulation');
InputDirectory = '../Input/GeoStat_ADM';
InputFile = strcat(InputDirectory, '/GeoStat.txt');
ResSimulator;
% disp('Start Bilinear interpolation')
% InputDirectory = '../Input/SPE10BRuns/Bilinear/DS03';
% InputFile = strcat(InputDirectory, '/SPE10B.txt');
% ResSimulator;
% clear all
% InputDirectory = '../Input/SPE10BRuns/Bilinear/DS02';
% InputFile = strcat(InputDirectory, '/SPE10B.txt');
% ResSimulator;
% InputDirectory = '../Input/SPE10BRuns/Bilinear/DS01';
% InputFile = strcat(InputDirectory, '/SPE10B.txt');
% ResSimulator;
% InputDirectory = '../Input/SPE10BRuns/Bilinear/DS007';
% InputFile = strcat(InputDirectory, '/SPE10B.txt');
% ResSimulator;
% InputDirectory = '../Input/SPE10BRuns/Bilinear/DS02';
% InputFile = strcat(InputDirectory, '/SPE10B.txt');
% ResSimulator;
% InputDirectory = '../Input/SPE10BRuns/Bilinear/DS003';
% InputFile = strcat(InputDirectory, '/SPE10B.txt');
% ResSimulator;
% InputDirectory = '../Input/SPE10BRuns/Bilinear/DS001';
% InputFile = strcat(InputDirectory, '/SPE10B.txt');
% ResSimulator;
% disp('end of Bilinear interpolation');
% disp('');
% disp('end of Bilinear interpolation');

%disp('Start MS interpolation')
% InputDirectory = '../Input/SPE10BRuns/MS/DS03';
% InputFile = strcat(InputDirectory, '/SPE10B.txt');
% ResSimulator;
% clear all
% InputDirectory = '../Input/SPE10BRuns/MS/DS02';
% InputFile = strcat(InputDirectory, '/SPE10B.txt');
% ResSimulator;
% InputDirectory = '../Input/SPE10BRuns/MS/DS01';
% InputFile = strcat(InputDirectory, '/SPE10B.txt');
% ResSimulator;
% InputDirectory = '../Input/SPE10BRuns/MS/DS007';
% InputFile = strcat(InputDirectory, '/SPE10B.txt');
% ResSimulator;
% InputDirectory = '../Input/SPE10BRuns/MS/DS02';
% InputFile = strcat(InputDirectory, '/SPE10B.txt');
% ResSimulator;
% InputDirectory = '../Input/SPE10BRuns/MS/DS003';
% InputFile = strcat(InputDirectory, '/SPE10B.txt');
% ResSimulator;
% InputDirectory = '../Input/SPE10BRuns/MS/DS001';
% InputFile = strcat(InputDirectory, '/SPE10B.txt');
% ResSimulator;
%disp('end of MS interpolation');

%disp('Start Constant interpolation');
% InputDirectory = '../Input/SPE10BRuns/Constant/DS03';
% InputFile = strcat(InputDirectory, '/SPE10B.txt');
% ResSimulator;
% clear all
% InputDirectory = '../Input/SPE10BRuns/Constant/DS02';
% InputFile = strcat(InputDirectory, '/SPE10B.txt');
% ResSimulator;
% InputDirectory = '../Input/SPE10BRuns/Constant/DS01';
% InputFile = strcat(InputDirectory, '/SPE10B.txt');
% ResSimulator;
% InputDirectory = '../Input/SPE10BRuns/Constant/DS007';
% InputFile = strcat(InputDirectory, '/SPE10B.txt');
% ResSimulator;
%InputDirectory = '../Input/SPE10BRuns/Constant/DS02';
%InputFile = strcat(InputDirectory, '/SPE10B.txt');
%ResSimulator;
% InputDirectory = '../Input/SPE10BRuns/Constant/DS003';
% InputFile = strcat(InputDirectory, '/SPE10B.txt');
% ResSimulator;
% InputDirectory = '../Input/SPE10BRuns/Constant/DS001';
% InputFile = strcat(InputDirectory, '/SPE10B.txt');
% ResSimulator;
%disp('end of Constant interpolation');

disp('End multiple runs');