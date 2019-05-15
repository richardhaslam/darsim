% LTS FIM strategy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef LTS_RefBc_enforcer < handle
    properties
        
    end
    methods
       function ComputeBcRefZone(obj, Refinement)
            % If the neighbour of a fine cell is in the coarse zone I need to
            % compute the flux at n+1 and use it as a Neumann b.c.
            Cells = find(Refinement);
            for c = Cells
                for nc = Grid.Neighbours(c, :)
                    if Ref(nc)
                        % Flux for b.c.
                        obj.LTSBcEvaluator()
                    end
                end
            end            
        end
    end
end