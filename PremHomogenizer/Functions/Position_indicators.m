%*********************************************************************
%*                   Indicators OF POSITION FOR EACH MESH            *
%*********************************************************************
%
% This function computes the propierties of each mesh
%
% Micro_geo.posGrad  = Numeration of degrees of freedom asociated to the 
% gradient for each triangle
% Micro_geo.posPres  =  = Numeration of degrees of freedom asociated to the 
% pressure for each triangle
% Micro_geo.posFlujo = Numeration of degrees of freedom asociated to the 
% numerical flux for each edge (triangle);
% Micro_geo.edgesKnum= Numeration of the edges of each triangle;
% Micro_geo.areaK    = Area of each triangle; 
% Micro_geo.normals  = normal vector of each vector; 
% Micro_geo.size     = mesh size (max(max(edgelenght)));
%
%***------------------------------------
% Manuela Bastidas - 2017.

function geo = Position_indicators(geo,gr)

%% MESH INFO
element      = geo.element;
coordinate   = geo.coordinate;
nElement     = geo.nElement;
nodes2edge   = geo.nodes2edge;

%% Indicators
% Vectorial position
posGrad  = cell(nElement,1);

% Scalar position
posPres    = cell(nElement,1);

% Numerical Flux position
posFlujo  = cell(nElement,1);

% Sort edges in a triangle
edgesKnum = cell(nElement,1);

% Area info
areaK = cell(nElement,1);

%normal vectors
normals = cell(nElement,1);

% Mesh size
hsize = zeros(nElement,1);

for elemento = 1:nElement
    
    % Coord (x;y) triangle vertex
    coord = coordinate(element(elemento,:),:)';
    
    % Unknows position / numeration
    posPres{elemento,1}  = 3*elemento-2:3*elemento;
    if gr == 0
        posGrad{elemento,1}  = [2*elemento-1, 2*elemento];
    else
        posGrad{elemento,1}  = 6*elemento-5:6*elemento;
    end
    
    % Numeration of edges for each triangle 
    edgesKnum{elemento,1} = diag(nodes2edge(element(elemento,[2 3 1]),...
                                 element(elemento,[3 1 2])));
    posFlujo{elemento,1}  = [2*edgesKnum{elemento,1}'-1;...
                             2*edgesKnum{elemento,1}'];

    % Triangle area
    areaK{elemento,1}  = abs(det([coord; ones(1,3)])/2);
    
    % coord egdes : [1:2] inicial , [3:4] final
    p      = [coord(:,[2 3 1]);coord(:,[3 1 2])];
    % Length edge
    le     = [norm(p(1:2,1)-p(3:4,1));norm(p(1:2,2)-p(3:4,2));...
              norm(p(1:2,3)-p(3:4,3))];
          
    for n=1:3
        % EN CONTRA DEL RELOJ -- ASI DEBE SER SIEMPRE
        % normalvector(1:2,n) = (p(1:2,n)-p(3:4,n))'*[0,1;-1,0]/le(n);
        normals{elemento,1}(1:2,n) = (p(1:2,n)-p(3:4,n))'*[0,1;-1,0]/le(n);
    end
    
    hsize(elemento) = max(le);

end

%% Micro_geo information
geo.posGrad  = posGrad;
geo.posPres  = posPres;
geo.posFlujo = posFlujo;
geo.edgesKnum= edgesKnum;
geo.areaK    = areaK; 
geo.normals  = normals; 
geo.size     = max(hsize);
