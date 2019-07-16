Case1.Dir = '../Input/ImmHomo_Fine/';
Case1.File = 'ImmHomo.txt';
Case2.Dir = '../Input/ImmHomo_LTS/';
Case2.File = 'ImmHomo.txt';
Error2VTKFile(Case1, Case2, '../', 10);

Fine = load(strcat(Case1.Dir, 'Output/SolverStats.txt'));
LTS = load(strcat(Case2.Dir, 'Output/LTSComplexity.txt'));

FineComplexity = sum(Fine(:, 3)) * 100*100;
LTSComplexity = sum(LTS);

disp(['Fine complexity: ', num2str(FineComplexity)]);
disp(['LTS complexity: ', num2str(LTSComplexity)]);