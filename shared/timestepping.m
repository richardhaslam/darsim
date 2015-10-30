function dt = timestepping(Fluid, S, Grid, U, Wells)
%Timestepping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Returns the timestep size based on the velocity field, the Grid and the
%fluid present. CFL condition is used. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nx=Grid.Nx;
Ny=Grid.Ny;
por=Grid.por;
pv=por*Grid.Volume;   %Void Volume in each cell
CFL=Grid.CFL;

%Compute wave speed
Sx=zeros(Nx+1,Ny);
Sy=zeros(Nx,Ny+1);
dfx=zeros(Nx+1,Ny);
dfy=zeros(Nx,Ny+1);
Sx(1,1) = (S(1,1) + 1)/2;
Sy(1,1) = (S(1,1) + 1)/2;
Sx(2:Nx,:)=(S(1:Nx-1,:)+S(2:Nx,:))/2;
Sy(:,2:Ny)=(S(:,1:Ny-1)+S(:,2:Ny))/2;
for i=1:Nx
    for j=1:Ny
        dfx(i,j)=Derivative(Sx(i,j),Fluid);
        dfy(i,j)=Derivative(Sy(i,j),Fluid);
    end
end
Lambdax=dfx.*U.x;
Lambday=dfy.*U.y;
Lambdax(1,1) = dfx(1,1)*max(max(Wells.Fluxes))/2;
Lambday(1,1) = dfy(1,1)*max(max(Wells.Fluxes))/2;
%Compute timestep size
dtx=CFL*pv/max(abs(Lambdax(:)));
dty=CFL*pv/max(abs(Lambday(:)));
dt=min(dtx,dty);

%Other way of computing the timestep size.
% XP=max(U.x,0); 
% XN=min(U.x,0);        % influx and outflux, x−faces
% YP=max(U.y,0); 
% YN=min(U.y,0);        % influx and outflux, y−faces
% fi=max(q,0);
% 
% Vi = XP(1:Nx,:)+YP(:,1:Ny)+...    % total flux into
%      -XN(2:Nx+1,:)-YN(:,2:Ny+1); % each gridblock
% 
% pm = min(pv./(Vi(:)+fi));            % estimate of influx
% dt = round(((1-Fluid.swc-Fluid.sor)/3)*pm); % CFL restriction
end