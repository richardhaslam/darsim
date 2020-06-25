% Table relative permeability model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 15 July 2016
%Last modified: 15 December 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef relperm_model_table < relperm_model
    properties
        Table
        Type
        Sw_from_Table
        kr_w_from_Table
        kr_nw_from_Table 
        krMax_w
        krMax_nw
        s_start
        s_end_w 
        s_end_nw
        n_w
        n_nw
        S_w
        S_nw
    end
    methods
        function obj = relperm_model_table(type,table)
            % Each table has 4 columns (Sw, Pc, Krw, Krnw).
            obj.Type = type;
            switch type
                case('Imbibition')
                    obj.Table.Data = table{1};
                    obj.S_irr(1) = min( obj.Table.Data(:,1) );
                    obj.S_irr(2) = 1 - max( obj.Table.Data(:,1) );
                case('Drainage')
                    obj.Table.Data = table{1};
                    obj.S_irr(1) = min( obj.Table.Data(:,1) );
                    obj.S_irr(2) = 1 - max( obj.Table.Data(:,1) );
                case('Cyclic')
                    %                     obj.Table.Imbibition.TableData = table{1};
                    %                     obj.Table.Draianage.TableData = table{2};
                    error('Cyclic rel perm is not implemented yet');
            end
            obj.PrepearingCoreyModel();
        end
        function PrepearingCoreyModel(obj)
            obj.Table.Data = sortrows(obj.Table.Data,1);
            obj.Sw_from_Table = obj.Table.Data(:,1);
            obj.kr_w_from_Table = obj.Table.Data(:,3);
            obj.kr_nw_from_Table = obj.Table.Data(:,4);
            
            obj.s_start = obj.Sw_from_Table(1);
            obj.s_end_w = obj.Sw_from_Table(end);
            % Rescale saturations of wetting
            S_w = (obj.Sw_from_Table-obj.s_start)/(1-obj.s_start-(1-obj.s_end_w));
            S_w = max(S_w, 0);
            obj.krMax_w = max(obj.kr_w_from_Table); 
            n_w_guess = 2;
            [obj.n_w, Rsquared_Corey_w] = fminsearch(PowerLawRelationship(n_w_guess,S_w,obj.kr_w_from_Table),);
               
            nonzeroKrnw = length(nonzeros(obj.kr_nw_from_Table)); 
            if nonzeroKrnw+2 <= length(obj.kr_nw_from_Table)
                obj.s_end_nw = obj.Sw_from_Table(nonzeroKrnw+2);
                
                % Rescale saturations of nonwetting
                S_nw = (obj.Sw_from_Table-obj.s_start)/(1-obj.s_start-(1-obj.s_end_nw));
                S_nw = max(S_nw, 0);
                S_nw = 1 - S_nw;
                obj.krMax_nw = max(obj.kr_nw_from_Table);
                n_nw_guess = 2;
                [obj.n_nw, Rsquared_Corey_nw] = fminsearch('PowerLawRelationship',n_nw_guess,[],S_nw(1:nonzeroKrnw+2),obj.kr_nw_from_Table(1:nonzeroKrnw+2));
            else
                obj.s_end_nw = obj.Sw_from_Table(nonzeroKrnw+1);
                
                % Rescale saturations of nonwetting
                S_nw = (obj.Sw_from_Table-obj.s_start)/(1-obj.s_start-(1-obj.s_end_nw));
                S_nw = max(S_nw, 0);
                S_nw = 1 - S_nw;
                obj.krMax_nw = max(obj.kr_nw_from_Table);
                n_nw_guess = 2;
                [obj.n_nw, Rsquared_Corey_nw] = fminsearch('PowerLawRelationship',n_nw_guess,[],S_nw(1:nonzeroKrnw+1),obj.kr_nw_from_Table(1:nonzeroKrnw+1));
            end
            
            % Phase 1 relative permeability
            kr(:,1) = obj.krMax_w * S_w.^obj.n_w;
            % Phase 2 relative permeability
            kr(:,2) = obj.krMax_nw * S_nw.^obj.n_nw;
            
            kr(obj.Sw_from_Table < obj.s_start, 1) = 0;
            kr(obj.Sw_from_Table < obj.s_start, 2) = 1;
            kr(obj.Sw_from_Table > obj.s_end_w, 1) = 1;
            kr(obj.Sw_from_Table > obj.s_end_nw, 2) = 0;
            
%             % Plot
%             figure
%             plot(obj.Sw_from_Table, obj.kr_w_from_Table, 'k*')
%             hold on
%             plot(obj.Sw_from_Table, obj.kr_nw_from_Table, 'k*')
%             hold on
%             plot(obj.Sw_from_Table, kr(:,1), 'k-')
%             hold on
%             plot(obj.Sw_from_Table, kr(:,2), 'k-o')
%             legend1 = strcat('krBrineCorey n= ', num2str(obj.n_w), '  r-square = ', num2str(1-Rsquared_Corey_w));
%             legend2 = strcat('krH2Corey n=', num2str(obj.n_nw), '  r-square = ', num2str(1-Rsquared_Corey_nw));
%             legend('SamplePointskrBrine', 'SamplePointskrH2',legend1,legend2);
%             legend('Location','west')
%             xlabel('Sw')
%             ylabel('kr')
%             xlim([0 1])
%             ylim([0 1])
%             title('RelPerm')
        end
        function kr = ComputeRelPerm(obj, ~, s)
            
            % Rescale saturations of wetting
            obj.S_w = (s-obj.s_start)/(1-obj.s_start-(1-obj.s_end_w));
            obj.S_w = max(obj.S_w, 0);
                
            % Rescale saturations of nonwetting
            obj.S_nw = (s-obj.s_start)/(1-obj.s_start-(1-obj.s_end_nw));
            obj.S_nw = max(obj.S_nw, 0);
            obj.S_nw = 1 - obj.S_nw; 
            
            % Phase 1 relative permeability
            kr(:,1) = obj.krMax_w * obj.S_w.^obj.n_w;
            % Phase 2 relative permeability
            kr(:,2) = obj.krMax_nw * obj.S_nw.^obj.n_nw;
            
            kr(s < obj.s_start, 1) = 0;
            kr(s < obj.s_start, 2) = 1;
            kr(s > obj.s_end_w, 1) = 1;
            kr(s > obj.s_end_nw, 2) = 0;
             
        end
        
        function dkr = ComputeDerivative(obj, ~, s)
            
            % Phase 1 relative permeability
            dkr(:,1) = obj.krMax_w / (1-obj.s_start-(1-obj.s_end_w)) * obj.n_w * obj.S_w.^(obj.n_w-1);
            % Phase 2 relative permeability
            dkr(:,2) = obj.krMax_nw / (-1*(1-obj.s_start-(1-obj.s_end_nw))) * obj.n_nw * (1-obj.S_nw).^(obj.n_nw-1);
            
            dkr(s < obj.s_start, 1) = 0;
            dkr(s < obj.s_start, 2) = 0;
            dkr(s > obj.s_end_w, 1) = 0;
            dkr(s > obj.s_end_nw, 2) = 0;
        end
        function ddkr = ComputeSecondDerivative(obj, ~, s)
            
            % Phase 1 relative permeability
            ddkr(:,1) = obj.krMax_w * (1-obj.s_start-(1-obj.s_end_w))^(-2) * obj.n_w * (obj.n_w -1) * obj.S_w.^(obj.n_w-2);
            % Phase 2 relative permeability
            ddkr(:,2) = obj.krMax_nw * (1-obj.s_start-(1-obj.s_end_nw))^(-2) * obj.n_nw * (obj.n_nw -1) * (1-obj.S_nw).^(obj.n_nw-2);
            
            ddkr(s < obj.s_start, 1) = 0;
            ddkr(s < obj.s_start, 2) = 0;
            ddkr(s > obj.s_end_w, 1) = 0;
            ddkr(s > obj.s_end_nw, 2) = 0;
            
        end
        function E = PowerLawRelationship(n,s,kr)
            fit = max(kr) * s.^n;
            
            Rsquared = (  length(s) * sum((fit.*kr)) - sum(fit) * sum(kr)  ) / ...
                       (  (length(s) * sum(fit.^2) - (sum(fit))^2 ) * ( length(s) * sum(kr.^2) - (sum(kr))^2 )  )^0.5;
                   
            E = 1 - Rsquared^2;
        end
    end
end