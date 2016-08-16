%  ADM operators handler base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 16 August 2016
%Last modified: 16 August 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef operators_handler < handle
    properties
        R
        Pp
        ADMRest
        ADMProlp
        ADMProls
    end
    methods
        function obj = operators_handler(n)
            obj.R = cell(1, n);
            obj.Pp = cell(1, n);
        end
        function MsR = MsRestriction(obj, FineGrid, CoarseGrid, Nf, Nc, level)
            %% MSFV Restriction Operator
            MsR = zeros(Nc, Nf);
            if (level == 1)
                for r=1:Nc
                    %coordinates of fine cells contained in the coarse block
                    Imin = CoarseGrid.I(r) - floor((CoarseGrid.CoarseFactor(1) - 1)/2);
                    Imax = CoarseGrid.I(r) + ceil((CoarseGrid.CoarseFactor(1) - 1)/2);
                    Jmin = CoarseGrid.J(r) - floor((CoarseGrid.CoarseFactor(2) - 1)/2);
                    Jmax = CoarseGrid.J(r) + ceil((CoarseGrid.CoarseFactor(2) - 1)/2);
                    i=Imin:Imax;
                    j=Jmin:Jmax;
                    [p,q] = meshgrid(i, j);
                    pairs = [p(:), q(:)];
                    %indexes of the fine cells
                    c = pairs(:,1) + (pairs(:,2)-1)*FineGrid.Nx;
                    %I make 1 those columns
                    MsR(r,c) = 1;
                end
            else
                for r=1:Nc
                    for c=1:Nf
                        if (r == FineGrid.Father(c, level))
                            %I make 1 those columns
                            MsR(r,c) = 1;
                        end
                    end
                end
            end
            MsR = sparse(MsR);
        end
    end
    methods (Abstract)
        obj = BuildStaticOperators(obj);
        obj = BuildADMOperators(obj);
    end
end
