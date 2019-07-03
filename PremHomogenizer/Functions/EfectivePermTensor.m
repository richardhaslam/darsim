%*********************************************************************
%*                                                                   *
%*              EFECTIVE Diffusion tensor A _ Upscaled               *
%*                                                                   *
%*********************************************************************
%
% This functions compute the efective difussion tensor.
%
%***------------------------------------
%***Inputs: Solution of the micro cell problems.


function [K_Efective] = EfectivePermTensor(Micro_geo,xx,yy,K_micro,Vel1,Vel2)

element      = Micro_geo.element;
coordinate   = Micro_geo.coordinate;

K_Efective = zeros(2,2);

for j = 1:Micro_geo.nElement
    % Coord (x;y) triangle vertex
    coord = coordinate(element(j,:),:)';
    
    % coord egdes : [1:2] inicial , [3:4] final
    p      = [coord(:,[2 3 1]);coord(:,[3 1 2])];
    % Length edge
    le     = [norm(p(1:2,1)-p(3:4,1)) norm(p(1:2,2)-p(3:4,2))...
        norm(p(1:2,3)-p(3:4,3))];
    
    I = diag(Micro_geo.nodes2edge(Micro_geo.element(j,[2 3 1]),...
        Micro_geo.element(j,[3 1 2])));
    signum = ones(1,3);
    signum((j==Micro_geo.edge2element(I,4)))= -1;
    
    bari = sum(coord,2)/3;
    %     area = det([1,1,1;coord])/2;
    
    %% Calculation of the velocity in the baricenter
    % See Raviart-Thomas formulation
    N=bari(:)*ones(1,3)-coord;
    P1=N*1/2*diag(signum)*diag(le)*Vel1(I);
    P2=N*1/2*diag(signum)*diag(le)*Vel2(I);
    
    aux1 = interp2(xx,yy,K_micro,bari(1),bari(2),'nearest');
    aux2 = interp2(xx,yy,K_micro,bari(1),bari(2),'nearest');
    K_Efective = K_Efective + ([aux1 0;0 aux2]*Micro_geo.area(j) +[P1 P2]);
    
end
K_Efective =  K_Efective./sum(Micro_geo.area);