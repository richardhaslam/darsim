%%%%Multiple Runs%%%%
disp('Start multiple runs');
% %clear all
% N = 1;
% Problem = 'SPE10T';
% Directories(1:2,:) = [strcat('../ResultsAndImages/PressureConstrained2/',Problem,'/ADM/DS03/MSbasisF_2Levels/');...
%     strcat('../ResultsAndImages/PressureConstrained2/',Problem,'/ADM/DS03/Bilinear_2Levels/')];
% LevelsMaxima = [2,2,3,3];
% Interpolators = {'MS'; 'Homogeneous';'Constant'; 'Homogeneous'};
% tolerance = 0.3;
% for i=1:N
%     Directory = Directories(i, :);
%     levelMax = LevelsMaxima(i);
%     InterpolationOption = Interpolators{i};
%     ResSimulator;
% end
% disp('End DS03');
clear all
N = 3;
Problem = 'SPE10T';
Directories(1:3,:) = [strcat('../ResultsAndImages/PressureConstrained2/',Problem,'/ADM/DS02/Constant_2Levels/');...
    strcat('../ResultsAndImages/PressureConstrained2/',Problem,'/ADM/DS02/Bilinear_2Levels/');...
    strcat('../ResultsAndImages/PressureConstrained2/',Problem,'/ADM/DS02/MSbasisF_2Levels/')];
LevelsMaxima = [2,2,2,3];
Interpolators = {'Constant'; 'Homogeneous';'MS'; 'Homogeneous'};
tolerance = 0.2;
for i=1:N
    Directory = Directories(i, :);
    levelMax = LevelsMaxima(i);
    InterpolationOption = Interpolators{i};
    ResSimulator;
end
disp('End DS02');
clear all
N = 1;
Problem = 'SPE10T';
Directories(1:4,:) = [strcat('../ResultsAndImages/PressureConstrained2/',Problem,'/ADM/DS01/MSbasisF_2Levels/');...
    strcat('../ResultsAndImages/PressureConstrained2/',Problem,'/ADM/DS01/Bilinear_2Levels/');...
    strcat('../ResultsAndImages/PressureConstrained2/',Problem,'/ADM/DS01/Constant_3Levels/');...
    strcat('../ResultsAndImages/PressureConstrained2/',Problem,'/ADM/DS01/Bilinear_3Levels/')];
LevelsMaxima = [2,2,3,3];
Interpolators = {'MS'; 'Homogeneous';'Constant'; 'Homogeneous'};
tolerance = 0.1;
for i=1:N
    Directory = Directories(i, :);
    levelMax = LevelsMaxima(i);
    InterpolationOption = Interpolators{i};
    ResSimulator;
end
disp('End DS01');
clear all
N = 1;
Problem = 'SPE10T';
Directories(1:4,:) = [strcat('../ResultsAndImages/PressureConstrained2/',Problem,'/ADM/DS007/MSbasisF_2Levels/');...
    strcat('../ResultsAndImages/PressureConstrained2/',Problem,'/ADM/DS007/Bilinear_2Levels/');...
    strcat('../ResultsAndImages/PressureConstrained2/',Problem,'/ADM/DS007/Constant_3Levels/');...
    strcat('../ResultsAndImages/PressureConstrained2/',Problem,'/ADM/DS007/Bilinear_3Levels/')];
LevelsMaxima = [2,2,3,3];
Interpolators = {'MS'; 'Homogeneous';'Constant'; 'Homogeneous'};
tolerance = 0.07;
for i=1:N
    Directory = Directories(i, :);
    levelMax = LevelsMaxima(i);
    InterpolationOption = Interpolators{i};
    ResSimulator;
end
disp('End DS007');
clear all
N = 1;
Problem = 'SPE10T';
Directories(1:4,:) = [strcat('../ResultsAndImages/PressureConstrained2/',Problem,'/ADM/DS005/MSbasisF_2Levels/');...
    strcat('../ResultsAndImages/PressureConstrained2/',Problem,'/ADM/DS005/Bilinear_2Levels/');...
    strcat('../ResultsAndImages/PressureConstrained2/',Problem,'/ADM/DS005/Constant_3Levels/');...
    strcat('../ResultsAndImages/PressureConstrained2/',Problem,'/ADM/DS005/Bilinear_3Levels/')];
LevelsMaxima = [2,2,3,3];
Interpolators = {'MS'; 'Homogeneous';'Constant'; 'Homogeneous'};
tolerance = 0.05;
for i=1:N
    Directory = Directories(i, :);
    levelMax = LevelsMaxima(i);
    InterpolationOption = Interpolators{i};
    ResSimulator;
end
disp('End DS005');
clear all
N = 1;
Problem = 'SPE10T';
Directories(1:4,:) = [strcat('../ResultsAndImages/PressureConstrained2/',Problem,'/ADM/DS003/MSbasisF_2Levels/');...
    strcat('../ResultsAndImages/PressureConstrained2/',Problem,'/ADM/DS003/Bilinear_2Levels/');...
    strcat('../ResultsAndImages/PressureConstrained2/',Problem,'/ADM/DS003/Constant_3Levels/');...
    strcat('../ResultsAndImages/PressureConstrained2/',Problem,'/ADM/DS003/Bilinear_3Levels/')];
LevelsMaxima = [2,2,3,3];
Interpolators = {'MS'; 'Homogeneous';'Constant'; 'Homogeneous'};
tolerance = 0.03;
for i=1:N
    Directory = Directories(i, :);
    levelMax = LevelsMaxima(i);
    InterpolationOption = Interpolators{i};
    ResSimulator;
end
disp('End DS003');
clear all
N = 1;
Problem = 'SPE10T';
Directories(1:4,:) = [strcat('../ResultsAndImages/PressureConstrained2/',Problem,'/ADM/DS001/MSbasisF_2Levels/');...
    strcat('../ResultsAndImages/PressureConstrained2/',Problem,'/ADM/DS001/Bilinear_2Levels/');...
    strcat('../ResultsAndImages/PressureConstrained2/',Problem,'/ADM/DS001/Constant_3Levels/');...
    strcat('../ResultsAndImages/PressureConstrained2/',Problem,'/ADM/DS001/Bilinear_3Levels/')];
LevelsMaxima = [2,2,3,3];
Interpolators = {'MS'; 'Homogeneous';'Constant'; 'Homogeneous'};
tolerance = 0.01;
for i=1:N
    Directory = Directories(i, :);
    levelMax = LevelsMaxima(i);
    InterpolationOption = Interpolators{i};
    ResSimulator;
end
disp('End multiple runs');