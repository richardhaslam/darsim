%%%%Multiple Runs%%%%
disp('Start multiple runs');
clear all
disp('Fine scale');
InputDirectory = '../Input/SPE10BRuns/FineScale/';
InputFile = strcat(InputDirectory, '/SPE10B.txt');
ResSimulator;

clear all
disp('Start Bilinear interpolation')
InputDirectory = '../Input/SPE10BRuns/Bilinear/DS03';
InputFile = strcat(InputDirectory, '/SPE10B.txt');
ResSimulator;
clear all
InputDirectory = '../Input/SPE10BRuns/Bilinear/DS02';
InputFile = strcat(InputDirectory, '/SPE10B.txt');
ResSimulator;
InputDirectory = '../Input/SPE10BRuns/Bilinear/DS01';
InputFile = strcat(InputDirectory, '/SPE10B.txt');
ResSimulator;
InputDirectory = '../Input/SPE10BRuns/Bilinear/DS007';
InputFile = strcat(InputDirectory, '/SPE10B.txt');
ResSimulator;
InputDirectory = '../Input/SPE10BRuns/Bilinear/DS005';
InputFile = strcat(InputDirectory, '/SPE10B.txt');
ResSimulator;
InputDirectory = '../Input/SPE10BRuns/Bilinear/DS003';
InputFile = strcat(InputDirectory, '/SPE10B.txt');
ResSimulator;
InputDirectory = '../Input/SPE10BRuns/Bilinear/DS001';
InputFile = strcat(InputDirectory, '/SPE10B.txt');
ResSimulator;
disp('end of Bilinear interpolation');
disp('');
disp('end of Bilinear interpolation');

disp('Start MS interpolation')
InputDirectory = '../Input/SPE10B/MS/DS03';
InputFile = strcat(InputDirectory, '/SPE10B.txt');
ResSimulator;
clear all
InputDirectory = '../Input/SPE10B/MS/DS02';
InputFile = strcat(InputDirectory, '/SPE10B.txt');
ResSimulator;
InputDirectory = '../Input/SPE10B/MS/DS01';
InputFile = strcat(InputDirectory, '/SPE10B.txt');
ResSimulator;
InputDirectory = '../Input/SPE10B/MS/DS007';
InputFile = strcat(InputDirectory, '/SPE10B.txt');
ResSimulator;
InputDirectory = '../Input/SPE10B/MS/DS005';
InputFile = strcat(InputDirectory, '/SPE10B.txt');
ResSimulator;
InputDirectory = '../Input/SPE10B/MS/DS003';
InputFile = strcat(InputDirectory, '/SPE10B.txt');
ResSimulator;
InputDirectory = '../Input/SPE10B/MS/DS001';
InputFile = strcat(InputDirectory, '/SPE10B.txt');
ResSimulator;
disp('end of MS interpolation');

disp('Start Constant interpolation');
InputDirectory = '../Input/SPE10BRuns/Constant/DS03';
InputFile = strcat(InputDirectory, '/SPE10B.txt');
ResSimulator;
clear all
InputDirectory = '../Input/SPE10BRuns/Constant/DS02';
InputFile = strcat(InputDirectory, '/SPE10B.txt');
ResSimulator;
InputDirectory = '../Input/SPE10BRuns/Constant/DS01';
InputFile = strcat(InputDirectory, '/SPE10B.txt');
ResSimulator;
InputDirectory = '../Input/SPE10BRuns/Constant/DS007';
InputFile = strcat(InputDirectory, '/SPE10B.txt');
ResSimulator;
InputDirectory = '../Input/SPE10BRuns/Constant/DS005';
InputFile = strcat(InputDirectory, '/SPE10B.txt');
ResSimulator;
InputDirectory = '../Input/SPE10BRuns/Constant/DS003';
InputFile = strcat(InputDirectory, '/SPE10B.txt');
ResSimulator;
InputDirectory = '../Input/SPE10BRuns/Constant/DS001';
InputFile = strcat(InputDirectory, '/SPE10B.txt');
ResSimulator;
disp('end of Constant interpolation');

disp('End multiple runs');