% Output writer FS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef output_writer_FS < output_writer
    properties
    end
    methods
        function obj = output_writer_FS(dir, problem, n_inj, n_prod, n_timers,  n_stats, n_comp)
            obj@output_writer(dir, problem, n_inj, n_prod, n_timers,  n_stats, n_comp);
        end
        function PlotSolution(obj, ProductionSystem, DiscretizationModel)
            obj.Plotter.PlotSolution(ProductionSystem, DiscretizationModel);
            obj.Plotter.PlotPermeability(ProductionSystem, DiscretizationModel);
            obj.Plotter.PlotPorosity(ProductionSystem, DiscretizationModel);
            obj.Plotter.VTKindex = obj.Plotter.VTKindex + 1; 
        end
        function WriteSummary(obj, Summary)
            obj.WriteWellsData(Summary.Time, Summary.WellsData, Summary.NumberTimeSteps);
            obj.WriteCouplingStats(Summary.CouplingStats, Summary.NumberTimeSteps);
        end
    end
end
       