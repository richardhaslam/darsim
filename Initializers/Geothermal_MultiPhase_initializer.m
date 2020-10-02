% InitializerSinglePhase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Geothermal_MultiPhase_initializer < initializer
    properties
    end
    methods
        function obj = Geothermal_MultiPhase_initializer(names, values)
            obj@initializer(names, values);
        end
        function ComputeInitialState(obj, ProductionSystem, FluidModel, Formulation, DiscretizationModel)
            disp('Started multiphase geothermal initialization (P-H)');
            
            if isnan(ProductionSystem.Reservoir.State.Properties('hTfluid').Value)
                P = ProductionSystem.Reservoir.State.Properties('P_2').Value;
                T = ProductionSystem.Reservoir.State.Properties('T').Value;
                H = ProductionSystem.Reservoir.State.Properties('hTfluid');
                H.Value = FluidModel.Phases(1).ComputeWaterEnthalpy(P, T);
                for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                    P = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties('P_2').Value;
                    T = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties('T').Value;
                    H = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties('hTfluid');
                    H.Value = FluidModel.Phases(1).ComputeWaterEnthalpy(P, T);
                end
            elseif isnan(ProductionSystem.Reservoir.State.Properties('T').Value)
                P = ProductionSystem.Reservoir.State.Properties('P_2').Value;
                H = ProductionSystem.Reservoir.State.Properties('hTfluid').Value;
                T = ProductionSystem.Reservoir.State.Properties('T');
                T.Value = FluidModel.Phases(1).ComputeWaterTemperature(P, H);
                for f = 1 : ProductionSystem.FracturesNetwork.NumOfFrac
                    P = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties('P_2').Value;
                    H = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties('hTfluid').Value;
                    T = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties('T');
                    T.Value = FluidModel.Phases(1).ComputeWaterTemperature(P, H);
                end
            else
                error('Both enthalpy and temperature cannot be unknown. Specify an input value for at least one of them.');
            end
            
            Formulation.ComputeProperties(ProductionSystem, FluidModel);
            
            % Output initial status:
            disp('Initial conditions:')
            
            % Initial status for reservoir:
            disp(['Pressure:         ', num2str(max(ProductionSystem.Reservoir.State.Properties('P_1').Value/1e5)), ' [bar]' ]);
            disp(['Enthalpy:         ', num2str(max(ProductionSystem.Reservoir.State.Properties('hTfluid').Value)), ' [J/Kg]']);
            disp(['Temperature:      ', num2str(max(ProductionSystem.Reservoir.State.Properties('T').Value      )), ' [K]'   ]);
            disp(['Water Saturation: ', num2str(max(ProductionSystem.Reservoir.State.Properties('S_1').Value    ))           ]);
            disp(['Steam Saturation: ', num2str(max(ProductionSystem.Reservoir.State.Properties('S_2').Value    ))           ]);
            disp('---------------------------------------------------------');
            
            disp(newline);
        end
    end
end