% Output writer for ADM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
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
            obj.Plotter.PlotSolution(ProductionSystem, DiscretizationModel);
            obj.Plotter.PlotPermeability(ProductionSystem, DiscretizationModel);
            obj.Plotter.PlotADMGrid(ProductionSystem, DiscretizationModel);
            if obj.basisfunctions
                obj.Plotter.PlotBasisFunctions(DiscretizationModel.FineGrid, ...
                    DiscretizationModel.CoarseGrid, DiscretizationModel.OperatorsHandler.ProlongationBuilders(1).P, DiscretizationModel.Nf, DiscretizationModel.Nc);
                obj.basisfunctions = false;
            end
            if obj.dynamicBF
                pressure = 0;
                saturation = 1;
                if (pressure)
                    % Pressure
                    obj.Plotter.PlotDynamicBasisFunctions(DiscretizationModel.ReservoirGrid, DiscretizationModel.OperatorsHandler.ADMProl{1})
                end
                if(saturation)
                    % Saturation
                    obj.Plotter.PlotSaturationInterpolator(DiscretizationModel.ReservoirGrid, DiscretizationModel.OperatorsHandler.ADMProl{2},...
                        DiscretizationModel.OperatorsHandler.ProlongationBuilders(2).Pdelta, DiscretizationModel.OperatorsHandler.ProlongationBuilders(2).Pdeltac);
                end
            end
            obj.Plotter.VTKindex = obj.Plotter.VTKindex + 1; 
        end
        function WriteSummary(obj, Summary)
            obj.WriteWellsData(Summary.Time, Summary.WellsData, Summary.NumberTimeSteps);
            obj.WriteCouplingStats(Summary.CouplingStats, Summary.NumberTimeSteps);
            obj.WriteADMStats(Summary.ADMStats, Summary.NumberTimeSteps);
        end
        function WriteADMStats(obj, ADMStats, Ndt)
            % Stats
            fileID = fopen(strcat(obj.Directory,'ADMStats.txt'),'w');
            [~, nc] = size(ADMStats);
            format = '%8d';
            for i =1:nc-2
                format = [format, ' %8d'];
            end
            format = [format, ' %8.2f\n'];
            fprintf(fileID, format, ADMStats(1:Ndt, :)');
            fclose(fileID);
        end
    end
end
