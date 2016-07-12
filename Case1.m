%%%%Case1 Runs%%%%
disp('Start multiple runs');
% 
Directories{1} = '../ResultsAndImages/ECMORPaper/Case1/MS/DS002';
Files{1} = '../ResultsAndImages/ECMORPaper/Case1/MS/DS002/GeoStat.txt';
Directories{2} = '../ResultsAndImages/ECMORPaper/Case1/MS/DS005';
Files{2} = '../ResultsAndImages/ECMORPaper/Case1/MS/DS005/GeoStat.txt';
Directories{3} = '../ResultsAndImages/ECMORPaper/Case1/MS/DS010';
Files{3} = '../ResultsAndImages/ECMORPaper/Case1/MS/DS010/GeoStat.txt';
Directories{4} = '../ResultsAndImages/ECMORPaper/Case1/MS/DS015';
Files{4} = '../ResultsAndImages/ECMORPaper/Case1/MS/DS015/GeoStat.txt';
Directories{5} = '../ResultsAndImages/ECMORPaper/Case1/MS/DS020';
Files{5} = '../ResultsAndImages/ECMORPaper/Case1/MS/DS020/GeoStat.txt';
for i=2:5
    InputDirectory = Directories{i};
    InputFile = Files{i};
    ResSimulator(Directories{i},Files{i});
end


Directories{1} = '../ResultsAndImages/ECMORPaper/Case1/Bilinear/DS002';
Files{1} = '../ResultsAndImages/ECMORPaper/Case1/Bilinear/DS002/GeoStat.txt';
Directories{2} = '../ResultsAndImages/ECMORPaper/Case1/Bilinear/DS005';
Files{2} = '../ResultsAndImages/ECMORPaper/Case1/Bilinear/DS005/GeoStat.txt';
Directories{3} = '../ResultsAndImages/ECMORPaper/Case1/Bilinear/DS010';
Files{3} = '../ResultsAndImages/ECMORPaper/Case1/Bilinear/DS010/GeoStat.txt';
Directories{4} = '../ResultsAndImages/ECMORPaper/Case1/Bilinear/DS015';
Files{4} = '../ResultsAndImages/ECMORPaper/Case1/Bilinear/DS015/GeoStat.txt';
Directories{5} = '../ResultsAndImages/ECMORPaper/Case1/Bilinear/DS020';
Files{5} = '../ResultsAndImages/ECMORPaper/Case1/Bilinear/DS020/GeoStat.txt';

for i=2:5
    InputDirectory = Directories{i};
    InputFile = Files{i};
    ResSimulator(Directories{i},Files{i});
end


Directories{1} = '../ResultsAndImages/ECMORPaper/Case1/Constant/DS002';
Files{1} = '../ResultsAndImages/ECMORPaper/Case1/Constant/DS002/GeoStat.txt';
Directories{2} = '../ResultsAndImages/ECMORPaper/Case1/Constant/DS005';
Files{2} = '../ResultsAndImages/ECMORPaper/Case1/Constant/DS005/GeoStat.txt';
Directories{3} = '../ResultsAndImages/ECMORPaper/Case1/Constant/DS010';
Files{3} = '../ResultsAndImages/ECMORPaper/Case1/Constant/DS010/GeoStat.txt';
Directories{4} = '../ResultsAndImages/ECMORPaper/Case1/Constant/DS015';
Files{4} = '../ResultsAndImages/ECMORPaper/Case1/Constant/DS015/GeoStat.txt';
Directories{5} = '../ResultsAndImages/ECMORPaper/Case1/Constant/DS020';
Files{5} = '../ResultsAndImages/ECMORPaper/Case1/Constant/DS020/GeoStat.txt';

for i=2:5
    InputDirectory = Directories{i};
    InputFile = Files{i};
    ResSimulator(Directories{i},Files{i});
end


disp('End multiple runs');