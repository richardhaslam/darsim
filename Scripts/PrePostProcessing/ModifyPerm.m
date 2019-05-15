%Modify Permeability files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Created: 2 May 2016
%Last modified: 5 Feb 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Angles = [0, 15, 30, 45, 98];
NewAngles = ['00deg'; '15deg'; '30deg'; '45deg'; '98deg'];

for i=1:5
    %Create new folder
    ReadingDirectory = strcat('../Permeability/',num2str(Angles(i)),' deg');
    WritingDirectory = strcat('../Permeability/', NewAngles(i,:));
    mkdir(WritingDirectory);
    for realization=1:20
        %Load it
        file = strcat(ReadingDirectory, '/Perm_', num2str(Angles(i)), '_', num2str(realization));
        data = load(file);
        %Modify it
        newvector = [data(1:2); 1; data(3:end)];
        %Save it in a new file
        Newfile = strcat(WritingDirectory, '/', 'Perm_',NewAngles(i,:),'_', num2str(realization), '.txt');
        disp(Newfile);
        save(Newfile, '-ascii', 'newvector');
    end
end
