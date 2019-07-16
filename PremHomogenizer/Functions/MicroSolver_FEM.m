%*********************************************************************
%*                                                                   *
%*                    Micro SOLVER                                   *
%*                                                                   *
%*********************************************************************
%
% This function assamble and solve the linear system associated with the
% micro-cell problems. 
% Output: Micro_geo: same + area 
%         Micro_Sol1&2: Gradient of the micro-scale solution

function [Micro_geo,Micro_sol1,Micro_sol2] = MicroSolver_FEM(Micro_geo,xx,yy,K_micro,Kx,Ky)

% Assemble matrices B and C
B = sparse(Micro_geo.noedges, Micro_geo.noedges);
C = sparse(Micro_geo.noedges,Micro_geo.nElement);

% Volume force
b1 = sparse(Micro_geo.noedges+size(Micro_geo.element,1),1);
b2 = sparse(Micro_geo.noedges+size(Micro_geo.element,1),1);

for j = 1:Micro_geo.nElement
    % Auxiliar variables
    coord = Micro_geo.coordinate(Micro_geo.element(j,:),:)';
    I = diag(Micro_geo.nodes2edge(Micro_geo.element(j,[2 3 1]),Micro_geo.element(j,[3 1 2])));
    signum = ones(1,3);
    signum((j==Micro_geo.edge2element(I,4)))= -1;
    
    bari = sum(coord,2)./3;
    aux1 = interp2(xx,yy,K_micro,bari(1),bari(2),'nearest');
    Micro_geo.area(j) = det([1,1,1;coord])/2;
    
    % Matrices of MFEM
    B(I,I)= B(I,I)+ diag(signum)*...
        stimaB(coord,[1/aux1 0;0 1/aux1])*diag(signum);
    
    n = coord(:,[3,1,2])-coord(:,[2,3,1]);
    C(I,j) = diag(signum)*[norm(n(:,1)) norm(n(:,2)) norm(n(:,3))]';
    
    % Sources of both micro-cell problems
    b1(Micro_geo.noedges+j)= det([1,1,1; Micro_geo.coordinate(Micro_geo.element(j,:),:)']) * ...
        interp2(xx,yy,Kx,bari(1),bari(2))/2;
    b2(Micro_geo.noedges+j)= det([1,1,1; Micro_geo.coordinate(Micro_geo.element(j,:),:)']) * ...
        interp2(xx,yy,Ky,bari(1),bari(2))/2;
end

% Global stiffness matrix
M = [B ,         C,      ;
    C', sparse(size(C,2),size(C,2))];

[x1,x2] = Micro_Boundary(M,b1,b2,Micro_geo);

Micro_sol1 = x1(1:end-Micro_geo.nElement);
Micro_sol2 = x2(1:end-Micro_geo.nElement);

end

function B=stimaB(coord,A)
N=coord(:)*ones(1,3)-repmat(coord,3,1);

D=diag([norm(N([5,6],2)) norm(N([1,2],3)) norm(N([1,2],2))]);

M=spdiags([ones(6,1),ones(6,1),2*ones(6,1),ones(6,1),ones(6,1)],...
    [-4,-2,0,2,4],6,6);

N_aux = (repmat(diag(A),3,3).*N)';
B = D*N_aux*M*N*D/(24*det([1,1,1;coord]));
end


