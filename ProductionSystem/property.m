classdef property < matlab.mixin.Copyable
    % This is a generic Property
    properties
        Value
        Valmin
        Valmax
        Plot
    end
    methods
        function obj = property(row, col, plot, minValue, maxValue, init)
            % Creates matrix with property values
            obj.Value = zeros(row, col);
            switch (nargin)
                case(3)
                    obj.Plot = plot;
                case(5)
                    obj.Valmin = minValue;
                    obj.Valmax = maxValue;
                case(6)
                    obj.Valmin = minValue;
                    obj.Valmax = maxValue;
                    obj.Value = init * ones(row, col);
            end
        end
        function update(obj, delta)
            obj.Value = obj.Value + delta;
        end
    end
end

