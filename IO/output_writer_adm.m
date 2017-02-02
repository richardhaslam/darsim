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
        basisfunctions
        dynamicBF
    end
    methods
        function obj = output_writer_adm(dir, problem, n_inj, n_prod, n_timers,  n_stats, n_comp)
            obj@output_writer(dir, problem, n_inj, n_prod, n_timers,  n_stats, n_comp);
            obj.basisfunctions = true;
            obj.dynamicBF = false;
        end
        function PlotSolution(obj, ProductionSystem, DiscretizationModel)
            obj.Plotter.PlotSolution(ProductionSystem.Reservoir.State, DiscretizationModel.ReservoirGrid);
            obj.Plotter.PlotPermeability(DiscretizationModel.ReservoirGrid, ProductionSystem.Reservoir.K);
            obj.Plotter.PlotADMGrid(DiscretizationModel.ReservoirGrid, DiscretizationModel.CoarseGrid);
            if obj.basisfunctions
                obj.Plotter.PlotBasisFunctions(DiscretizationModel.ReservoirGrid, DiscretizationModel.CoarseGrid, DiscretizationModel.OperatorsHandler.Pp);
                obj.basisfunctions = false;
            end
            if obj.dynamicBF
                obj.Plotter.PlotDynamicBasisFunctions(DiscretizationModel.ReservoirGrid, DiscretizationModel.OperatorsHandler.ADMProlp)
            end
        end
        function WriteSummary(obj, Summary)
            obj.WriteWellsData(Summary.Time, Summary.WellsData, Summary.NumberTimeSteps);
            obj.WriteCouplingStats(Summary.CouplingStats, Summary.NumberTimeSteps);
            obj.WriteADMStats(Summary.ADMStats, Summary.NumberTimeSteps);
        end
        function WriteADMStats(obj, ADMStats, Ndt)
             %Stats
            fileID = fopen(strcat(obj.Directory,'ADMStats.txt'),'w');
            fprintf(fileID, '%8d %8d %8d %8.2f\n', ADMStats(1:Ndt, :)');
            fclose(fileID);
        end
    end
end
