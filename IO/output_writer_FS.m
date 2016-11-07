% Output writer FS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 6 September 2016
%Last modified: 6 September 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef output_writer_FS < output_writer
    properties
    end
    methods
        function obj = output_writer_FS(dir, problem, n_inj, n_prod, n_timers,  n_stats)
            obj@output_writer(dir, problem, n_inj, n_prod, n_timers,  n_stats);
        end
        function PlotSolution(obj, ProductionSystem, DiscretizationModel)
            obj.Plotter.PlotSolution(ProductionSystem.Reservoir.State, DiscretizationModel.ReservoirGrid);
            obj.Plotter.PlotPermeability(DiscretizationModel.ReservoirGrid, ProductionSystem.Reservoir.K);
        end
        function WriteSummary(obj, Summary)
            obj.WriteWellsData(Summary.Time, Summary.WellsData, Summary.NumberTimeSteps);
            obj.WriteCouplingStats(Summary.CouplingStats, Summary.NumberTimeSteps);
        end
    end
end

   
       