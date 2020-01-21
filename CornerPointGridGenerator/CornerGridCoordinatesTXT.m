%% CornerPointGrid Converter for DARSim2
% Based on MRST, it reads the input files of CornerPointGrids (from Eclipse formating)
% and generates the input file that can be read by DARSim2.
% 
% Author: Janio Piguave
% Modified on: 2020/01/21
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
clc;

grdecl        = readGRDECL('NPD5.grdecl');
% Build Corner Point Nodes of each one of the cells
[x,y,z] = buildCornerPtNodes(grdecl);
Nx = grdecl.cartDims(1); Ny = grdecl.cartDims(2); Nz = grdecl.cartDims(3);
CellsN = linspace(1, Nx*Ny*Nz, Nx*Ny*Nz)';

% actnum = grdecl.ACTNUM;
% grdecl.ACTNUM = ones(prod(grdecl.cartDims),1);
% G = processGRDECL(grdecl, 'checkgrid', false);
G = processGRDECL(grdecl, 'Verbose', true);
G = computeGeometry(G);
plotGrid(G)
% plotGrid(G, ~actnum(G.cells.indexMap), 'FaceColor','none','EdgeAlpha',.1);
% plotGrid(G, actnum(G.cells.indexMap), 'FaceColor','y', 'EdgeAlpha',.1);
axis off; view(-155,80); zoom(1.7);

% plotGrid(G,'FaceColor','b','FaceAlpha', 0.1);
% 
% %%%%%%%%%%%%%
% % start: add add add
% c_cent = G.cells.centroids;
% f_cent = G.faces.centroids;
% coords = G.nodes.coords;
% 
% % Add circles around the centroids of each cell
% hold on;
% pltarg = {'MarkerSize',20,'LineWidth',2,'MarkerFaceColor',[.95 .95 .95]};
% 
% plot3(c_cent(:,1), c_cent(:,2), c_cent(:,3),'or',pltarg{:});
% 
% % % Plot triangles around face centroids
% % plot3(f_cent(:,1), f_cent(:,2), f_cent(:,3),'sg',pltarg{:});
% % 
% % % Plot squares around nodes
% % plot3(coords(:,1), coords(:,2), coords(:,3),'db',pltarg{:});
% 
% legend({'Grid', 'Cell', 'Face', 'Node'}, 'Location', 'SouthOutside', 'Orientation', 'horizontal')
% 
% % Plot cell/face centroids and nodes
% txtargs = {'FontSize',12,'HorizontalAlignment','left'};
% text(c_cent(:,1)-0.005, c_cent(:,2), c_cent(:,3), num2str((1:G.cells.num)'),txtargs{:});
% % text(f_cent(:,1)-0.0105, f_cent(:,2), f_cent(:,3), num2str((1:G.faces.num)'),txtargs{:});
% % text(coords(:,1)-0.012, coords(:,2), coords(:,3), num2str((1:G.nodes.num)'),txtargs{:});
% 
% title('Grid structure')
% hold off;
% 
% % end add add
% %%%%%%%%%%%%%%%%

R = size(x);
X = reshape(permute(reshape(x,R(1),2,[],R(3)),[1,3,2,4]),[],2,R(3));

A=[]; a=1;
for k=2:2:size(X,3)
    Aa(:,:,a)=[X(:,:,k-1),X(:,:,k)];                              %horizontal cat
    a=a+1;
end
for k=1:size(X,3)/2
    A=[A;Aa(:,:,k)];                                                 %vertical cat
end
A=A';
A=reshape(A,8,[])';

S = size(y);
Y = reshape(permute(reshape(y,S(1),2,[],S(3)),[1,3,2,4]),[],2,S(3));
B=[]; b=1;
for k=2:2:size(Y,3)
    Bb(:,:,b)=[Y(:,:,k-1),Y(:,:,k)];                              %horizontal cat
    b=b+1;
end
for k=1:size(Y,3)/2
    B=[B;Bb(:,:,k)];                                                 %vertical cat
end
B=B';
B=reshape(B,8,[])';

T = size(z);
Z = reshape(permute(reshape(z,T(1),2,[],T(3)),[1,3,2,4]),[],2,T(3));
C=[]; c=1;
for k=2:2:size(Z,3)
    Cc(:,:,c)=[Z(:,:,k-1),Z(:,:,k)];                              %horizontal cat
    c=c+1;
end
for k=1:size(Z,3)/2
    C=[C;Cc(:,:,k)];                                                 %vertical cat
end
C=C';
C=reshape(C,8,[])';

% A = x, B = y, C = z, 
% Matrix with the Nodes for the Vectors: SW_B, SE_B, NW_B, SW_T
% vector a = SE_B - SW_B
va = [(A(:,7) - A(:,3)), (B(:,7) - B(:,3)), (C(:,7) - C(:,3))];
% vector b = NW_B - SW_B
vb = [(A(:,4) - A(:,3)), (B(:,4) - B(:,3)), (C(:,4) - C(:,3))];
% vector c = NW_B - SW_B
vc = [(A(:,3) - A(:,1)), (B(:,3) - B(:,1)), (C(:,3) - C(:,1))];

f2cn = gridCellNo(G);

% Section 1
% Nodes_XYZ + Centroids + Vectors a b c
% Only Active Cells (Based on ACTNUM info)
Nodes_XYZ = [CellsN, grdecl.ACTNUM, A(:,2), B(:,2), C(:,2), A(:,6), B(:,6), C(:,6), A(:,1), B(:,1), C(:,1), A(:,5), B(:,5), C(:,5),...
    A(:,4), B(:,4), C(:,4), A(:,8), B(:,8), C(:,8), A(:,3), B(:,3), C(:,3), A(:,7), B(:,7), C(:,7), va, vb, vc];
TF = Nodes_XYZ(:,2) == 0;
Nodes_XYZ(TF,:) = []; Nodes_XYZ(:,2) = [];
Nodes_XYZi = [Nodes_XYZ(:,1:31), G.cells.centroids, Nodes_XYZ(:,32:34)];

% Nodes_XYZi = [A(:,2), B(:,2), C(:,2), A(:,6), B(:,6), C(:,6), A(:,1), B(:,1), C(:,1), A(:,5), B(:,5), C(:,5),...
%     A(:,4), B(:,4), C(:,4), A(:,8), B(:,8), C(:,8), A(:,3), B(:,3), C(:,3), A(:,7), B(:,7), C(:,7),...
%     G.cells.centroids, va, vb, vc];

% Section 2: CELL DATA
% dx - dy - dz
% dx = SE_B_x - SW_B_x  |  dy = NW_B_y - SW_B_y  |  dz = SW_B_z - SW_T_z
delta = [((A(:,7)-A(:,3))), ((B(:,4)-B(:,3))), ((C(:,3)-C(:,1)))];
Section2a = [CellsN, grdecl.ACTNUM, delta];
TF = Section2a(:,2) == 0;
Section2a(TF,:) = []; Section2a(:,2) = [];
Section2b = [Section2a, G.cells.volumes];

% Section 3: Faces Connected
% %%%% For getting the Array - Cells that belong a to a i Face
% % Cells N times (later related with faces & tags)
% f2cn = gridCellNo(G);
% % Cells + Faces + Tags
% F = [f2cn G.cells.faces];
% % Cell Index links to Faces Index 
% CtoF = accumarray(F(:,2), F(:,1), [], @(x){x});
% % make sure all elements are column vectors 
% % CtoF2 = cellfun(@(x) reshape(x,[],1), CtoF,'un',0); 
% CtoF2 = padcat(CtoF{:})';
% %%%%%%%%%%%%%%

% Create Face Index Vector (Total Number of Faces)
FacesN = linspace(1, G.faces.num, G.faces.num)';
% Assembly Matrix with all the Faces
% FC = [FacesN, G.faces.neighbors, CtoF2, G.faces.areas];
FC = [FacesN, G.faces.neighbors, G.faces.areas];
TF = (FC(:,2) == 0)|(FC(:,3) == 0);
FC(TF,:) = [];

% Section 4: Faces Disconnected
FD = [FacesN, G.faces.neighbors];
TF = (FD(:,2) ~= 0)&(FD(:,3) ~= 0);
FD(TF,:) = [];

%% Writing the output file
Directory = 'C:\Users\Janio Paul\DARSim2\MSRT\mrst-2019a_zip\'; 
OutputFileName = 'Test_3.txt';
disp( '******************* Writing the data into output text file *********************' );
% disp(['Deleting file ', OutputFileName]);
% delete( strcat(Directory,OutputFileName) );
disp(['Writing into file ', OutputFileName]);
fid = fopen(strcat(Directory,OutputFileName) , 'w+' );

fprintf(fid, '%% Node Locations for each one of the cells, Nx * Ny *Nz\n');
fprintf(fid, '%% Coordinates X,Y,Z of the eight nones for each one of the cells\n');
fprintf(fid, '\n');
fprintf(fid, '%% The grid resolution of reservoir:\n');
fprintf(fid, 'RESERVOIR_GRID   %8.0f x %8.0f x %8.0f [ - ]\n' , grdecl.cartDims(1) , grdecl.cartDims(2) , grdecl.cartDims(3) );
fprintf(fid, '\n');
fprintf(fid, '%% Coordinates of grid points arranged row-wise (through length) then column-wise (through columns)\n');
fprintf(fid,'%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %7s %7s %7s %7s %7s %7s %8s %8s %8s %7s %7s %7s\n', ...
          'C  ', ... 
          'NW_T(x)', 'NW_T(y)', 'NW_T(z)', ' NE_T(x)', 'NE_T(y)', 'NE_T(z)', ' SW_T(x)', 'SW_T(y)', 'SW_T(z)', ' SE_T(x)', 'SE_T(y)', 'SE_T(z)', ...
          ' NW_B(x)', 'NW_B(y)', 'NW_B(z)', ' NE_B(x)', 'NE_B(y)', 'NE_B(z)', ' SW_B(x)', 'SW_B(y)', 'SW_B(z)', ' SE_B(x)', 'SE_B(y)', 'SE_B(z)', ...
          ' CC(x)', 'CC(y)', 'CC(z)', ' a(x)', 'a(y)', 'a(z)', ' b(x)', 'b(y)', 'b(z)', ' c(x)', 'c(y)', 'c(z)');      
      
for ii = 1:size(Nodes_XYZi,1)
    fprintf(fid,'%u: % .5f %.5f %.5f % .5f %.5f %.5f % .5f %.5f %.5f % .5f %.5f %.5f % .5f %.5f %.5f % .5f %.5f %.5f % .5f %.5f %.5f % .5f %.5f %.5f % .5f %.5f %.5f % .5f %.5f %.5f % .5f %.5f %.5f % .5f %.5f %.5f\n', Nodes_XYZi(ii,:)');
end
% for ii = 1:size(Nodes_XYZ,1)
%  fprintf(fid, '%u: ', CellsN(ii));
%     fprintf(fid,'% .5f %.5f %.5f % .5f %.5f %.5f % .5f %.5f %.5f % .5f %.5f %.5f % .5f %.5f %.5f % .5f %.5f %.5f % .5f %.5f %.5f % .5f %.5f %.5f % .5f %.5f %.5f % .5f %.5f %.5f % .5f %.5f %.5f % .5f %.5f %.5f\n', Nodes_XYZ(ii,:)');
% end
fprintf(fid, '\n\n');
fprintf(fid, '%% Section 2: Cells Information\n');
fprintf(fid, '%% [Cell Index] + [dx] + [dy] + [dz] + [Cell Volume]:\n');
fprintf(fid, '\n');
fprintf(fid,'%s %5s %7s %7s %9s\n','C', 'dx', 'dy', 'dz', 'Volume');      

for ii = 1:size(Section2b,1)
    fprintf(fid,'%u %.5f %.5f %.5f %.5f\n', Section2b(ii,:)');
end

fprintf(fid, '\n\n');
fprintf(fid, '%% Section 3: Faces Connected Information\n');
fprintf(fid, '%% [Face Index] + [Cell Index] + [Cell Index] + [Cell Area]:\n');
fprintf(fid, '\n');
fprintf(fid,'%3s %3s %3s %5s\n','F','C1','C2','Area');      

for ii = 1:size(FC,1)
    fprintf(fid,'%3u %3u %3u %.4f\n', FC(ii,:)');
end

fprintf(fid, '\n\n');
fprintf(fid, '%% Section 4: Faces Disconnected Information\n');
fprintf(fid, '%% [Face Index] + [Cell Index] + [Cell Index]:\n');
fprintf(fid, '\n');
fprintf(fid,'%3s %3s %3s\n','F','C1','C2');      

for ii = 1:size(FD,1)
    fprintf(fid,'%3u %3u %3u\n', FD(ii,:)');
end