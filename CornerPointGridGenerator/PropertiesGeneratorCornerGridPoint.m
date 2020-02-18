close all; clear; clc;
% Read subset of ECLIPSE GRID file. Return: CartDims + ZCORN + COORD + ACTNUM
grdecl = readGRDECL('NPD5.grdecl');
G = processGRDECL(grdecl, 'Verbose', true);
% Compute geometry information (centroids, volumes, areas) of the cells
G = computeGeometry(G);
% Vector with the total number of cells
NumberCells = linspace(1, grdecl.cartDims(1)*grdecl.cartDims(2)*grdecl.cartDims(3), grdecl.cartDims(1)*grdecl.cartDims(2)*grdecl.cartDims(3))';
% WHY DOUBLE
ActiveCells = double([NumberCells grdecl.ACTNUM]);
LI = ActiveCells(:,2) == 0;                            % Logical Index                 
ActiveCells(LI,:) = [];                                % Delete Cells that are not active
ActiveCells(:,2) = [];                                 % Delete Columns of Active Cells
% Load Inputs Files
p = load('NPD5_Porosity.txt')';
K = load('NPD5_Permeability.txt')';
% Just Considered Active Cells
p = p(G.cells.indexMap);
poro = p;
poro_text = [ActiveCells poro];
K = K(G.cells.indexMap);
% Convert K values a diferent units
K = K .* milli * darcy;
perm = bsxfun(@times, [1 1 0.1], K);
perm_txt = [ActiveCells perm];
rock = makeRock(G, perm, poro);
[K, i, j] = permTensor(rock, G.griddim);
K_text = [ActiveCells K];

%% OUTPUT FILE: WRITING
Directory = 'C:\Users\Janio Paul\DARSim2\MSRT\mrst-2019a_zip\'; 
OutputFileName = 'CornerPointGrid_DARSim_RockPropertiesData.txt';
disp( '******************* Writing the data into output text file *********************' );
disp(['Writing into file ', OutputFileName]);

fid = fopen(strcat(Directory, OutputFileName), 'w+');
fprintf(fid, '%% Rock Properties of the Johansen Formation\n');
fprintf(fid, '%% Rock Properties Values of Active Cells of the Corner Point Grid Model\n');
fprintf(fid, '%% Permeability tensor is assumed to be diagonal\n');
fprintf(fid, '%% Vertical permeability (Kz) equaling one-tenth of the horizontal permeability (Kx, Ky)\n');
fprintf(fid, '%% Only the x-component Kx is given in the data file. Kx = Ky. Kz = 0.1 Kx\n');
fprintf(fid, '\n');
fprintf(fid, '%% The Grid Resolution of the Reservoir:\n');
fprintf(fid, 'RESERVOIR_GRID_NX\n');
fprintf(fid, '%d\n', grdecl.cartDims(1));
fprintf(fid, 'RESERVOIR_GRID_NY\n');
fprintf(fid, '%d\n', grdecl.cartDims(2));
fprintf(fid, 'RESERVOIR_GRID_NZ\n');
fprintf(fid, '%d\n', grdecl.cartDims(3));
fprintf(fid, 'ACTIVE_CELLS\n');
fprintf(fid, '%d\n', G.cells.num);

fprintf(fid, '\n');
fprintf(fid, '%% Section 1: Porosity - [Cell Index] + [Porosity]\n');
fprintf(fid, '%% Units: Fraction\n');
fprintf(fid, 'POROSITY\n');
 
fprintf(fid,'%s %s\n', 'Cell No.  ','Porosity');   
      
for ii = 1:size(poro_text,1)
    fprintf(fid,'%6.0d   ,  %.4f\n', poro_text(ii,:)');
end

fprintf(fid, '\n');
fprintf(fid, '%% Section 2: Pemerbility - [Cell Index] + [Permeability (Kx, Ky, Kz)]\n');
fprintf(fid, '%% Units: Metric\n');
fprintf(fid, 'PERMEABILITY_XYZ\n');
 
fprintf(fid,'%s %7s %10s %10s\n', 'Cell No.  ','Kx','Ky','Kz');   
   
for ii = 1:size(perm_txt,1)
    fprintf(fid,'%6.0d   ,  %.4e;%.4e;%.4e\n', perm_txt(ii,:)');
end

fprintf(fid, '\n');
fprintf(fid, '%% Section 3: Pemerbility - [Cell Index] + [Permeability (Kx, Ky, Kz)]\n');
fprintf(fid, '%% Units: Metric\n');
fprintf(fid, 'PERMEABILITY_TENSOR\n');
 
fprintf(fid,'%s %8s %10s %10s %14s %10s %10s %14s %10s %10s\n', 'Cell No.  ','Kxx','Kxy','Kxz','Kyx','Kyy','Kyz','Kzx','Kzy','Kzz');   
   
for ii = 1:size(K_text,1)
    fprintf(fid,'%6.0d   ,  %.4e;%.4e;%.4e  ,  %.4e;%.4e;%.4e  ,  %.4e;%.4e;%.4e\n', K_text(ii,:)');
end