%  Basis functions updater
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 16 August 2016
%Last modified: 16 August 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef bf_updater_ms < bf_updater
    properties
    end
    methods
        function ConstructPressureSystem(obj, FineGrid)
            %Assemble pressure matrix
            N =  FineGrid.N;
            Nx = FineGrid.Nx;
            Ny = FineGrid.Ny;
            
            %Construct pressure matrix
            x1 = reshape(Tx(1:Nx,:),N,1);
            x2 = reshape(Tx(2:Nx+1,:),N,1);
            y1 = reshape(Ty(:,1:Ny),N,1);
            y2 = reshape(Ty(:,2:Ny+1),N,1);
            DiagVecs = [-y2,-x2,y2+x2+y1+x1,-x1,-y1];
            DiagIndx = [-Nx,-1,0,1,Nx];
            obj.A = spdiags(DiagVecs,DiagIndx,N,N);
        end
        function ComputeBasisFunctions(obj)
        end
    end
end