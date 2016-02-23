%Time-step selection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dt = timestepping(Fluid, S, Grid, U, Wells)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Returns the timestep size based on the velocity field, the Grid and the
%fluid present. CFL condition is used. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nx=Grid.Nx;
Ny=Grid.Ny;
por=Grid.por;
pv=por*Grid.Volume;   %Void Volume in each cell
CFL=Grid.CFL;

%I take the worst possible scenario
s = Fluid.swc:0.01:1-Fluid.sor;
df = Derivative (s, Fluid);
dfmax = max(df);
Uxmax = max(max(U.x));
Uymax = max(max(U.y));
Lambdax = dfmax * Uxmax;
Lambday = dfmax * Uymax;
%Compute timestep size
dtx=CFL*pv/Lambdax;
dty=CFL*pv/Lambday;
dt=min(dtx,dty);

% %%%%%%%%%%%%%% Foam? Compute based on worst-case scenario
% 
% if (strcmp(Fluid.RelPerm,'Foam') == 1) 
% 	A = linspace(Fluid.swc,1-Fluid.sor,Ny);
% 	S = repmat(A,Nx,1);
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %Compute wave speed
% Sx = zeros(Nx+1,Ny);
% Sy = zeros(Nx,Ny+1);
% dfx = zeros(Nx+1,Ny);
% dfy = zeros(Nx,Ny+1);
% Sx(1,1) = (S(1,1) + 1)/2;
% Sy(1,1) = (S(1,1) + 1)/2;
% Sx(2:Nx,:) = (S(1:Nx-1,:)+S(2:Nx,:))/2;
% Sy(:,2:Ny) = (S(:,1:Ny-1)+S(:,2:Ny))/2;
% 
% tstart = tic;
% for i=1:Nx
%     for j=1:Ny
%         dfx(i,j) = Derivative(Sx(i,j),Fluid);
%         dfy(i,j) = Derivative(Sy(i,j),Fluid);
%     end
% end
% timesteptimer = toc(tstart)
% Lambdax=dfx.*U.x;
% Lambday=dfy.*U.y;
% Lambdax(1,1) = dfx(1,1)*max(max(Wells.Fluxes))/2;
% Lambday(1,1) = dfy(1,1)*max(max(Wells.Fluxes))/2;
% %Compute timestep size
% dtx=CFL*pv/max(abs(Lambdax(:)));
% dty=CFL*pv/max(abs(Lambday(:)));
% dt=min(dtx,dty);

% tstart = tic;
% %Other way of computing the timestep size.
% XP=max(U.x,0); 
% XN=min(U.x,0);        % influx and outflux, x−faces
% YP=max(U.y,0); 
% YN=min(U.y,0);        % influx and outflux, y−faces
% q = reshape(Wells.Fluxes, Nx*Ny, 1);
% fi=max(q,0);
% 
% Vi = XP(1:Nx,:)+YP(:,1:Ny)+...    % total flux into
%      -XN(2:Nx+1,:)-YN(:,2:Ny+1); % each gridblock
% 
% pm = min(pv./(Vi(:)+fi));            % estimate of influx
% dt = round(((1-Fluid.swc-Fluid.sor)/3)*pm); % CFL restriction
% option = toc(tstart)
end