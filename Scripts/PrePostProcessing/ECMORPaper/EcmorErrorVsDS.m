% ECMOR PAPER: ERRORS vs DS
Directory = '../../ResultsAndImages/ECMORPaper/Case3/';
%% Read data from files
%Fine Scale
file = strcat(Directory,'FineScale/Output/FIMPressure.txt');
Pf = load(file);
file = strcat(Directory,'FineScale/Output/FIMSaturation.txt');
Sf = load(file);
%ADM with MS interpolation
DS{1} = 'DS005';
DS{2} = 'DS010';
DS{3} = 'DS015';
DS{4} = 'DS020';
for i=1:4
    file = strcat(Directory,'MS/',DS{i},'/Output/ADMPressure.txt');
    ADM{i}.P = load(file);
    file = strcat(Directory,'MS/',DS{i},'/Output/ADMSaturation.txt');
    ADM{i}.S = load(file);
end


%% Compute errors
for i=1:4
    PressureError = zeros(50,1);
    SaturationError = zeros(50,1);
    for j=1:50
    PressureError (j) = norm(ADM{i}.P(:,j) - Pf(:,j), 2) / norm(Pf(:,j),2);
    SaturationError(j) = norm(ADM{i}.S(:,j) - Sf(:,j),2)/9801;
    end
    ADM{i}.errorP =  mean(PressureError);
    ADM{i}.errorS = mean(SaturationError);
end


%% Print tables for Latex
ds = [0.05;0.1;0.15;0.20];
errorp = zeros(4,1);
errors = zeros(4,1);
for i =1:4
errorp(i) = ADM{i}.errorP;
errors(i) = ADM{i}.errorS;
end
fileID = fopen(strcat(Directory, 'Case3ErrorP.txt'), 'w');
Print = [ds,errorp];
fprintf(fileID, '%5.3f %10.6f\n', Print');
fclose(fileID);
Print = [ds,errors];
fileID = fopen(strcat(Directory, 'Case3ErrorS.txt'), 'w');
fprintf(fileID, '%5.3f %10.6f\n', Print');
fclose(fileID);