function C = inputdataDARSim(z)
T = size(z);
Z = reshape(permute(reshape(z,T(1),2,[],T(3)),[1,3,2,4]),[],2,T(3));
C=[]; c=1;
for k=2:2:size(Z,3)
    Cc(:,:,c)=[Z(:,:,k-1),Z(:,:,k)];                             
    c=c+1;
end
for k=1:size(Z,3)/2
    C=[C;Cc(:,:,k)];                                               
end
C=C';
C_VTK1 = reshape(C,1,[])';
C=reshape(C,8,[])';
end