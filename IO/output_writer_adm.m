% Output writer for ADM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 6 September 2016
%Last modified: 6 September 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef output_writer_adm < output_writer
    properties
        FormatADM
    end
    methods
        function obj = output_writer_adm(dir, problem, n_inj, n_prod, n_timers,  n_stats)
            obj@output_writer(dir, problem, n_inj, n_prod, n_timers,  n_stats);
        end
        function PlotSolution(obj, ProductionSystem, DiscretizationModel)
            obj.Plotter.PlotSolution(ProductionSystem.Reservoir.State, DiscretizationModel.ReservoirGrid);
            obj.Plotter.PlotADMGrid(DiscretizationModel.ReservoirGrid, DiscretizationModel.CoarseGrid);
        end
        function WriteSummary(obj, Summary)
            obj.WriteWellsData(Summary.Time, Summary.WellsData, Summary.NumberTimeSteps);
            obj.WriteCouplingStats(Summary.CouplingStats, Summary.NumberTimeSteps);
            obj.WriteADMStats();
        end
        function WriteADMStats(obj)
        end
    end
end
