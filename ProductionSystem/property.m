classdef property < matlab.mixin.Copyable
    % This is a generic Property
    properties
        Value
        Type
        Valmin
        Valmax
        Plot
    end
    methods
        function obj = property(row, col, type, plot, minValue, maxValue, init)
            % Creates matrix with property values
            obj.Value = zeros(row, col);
            obj.Type = type;
            switch (nargin)
                case(4)
                    obj.Plot = plot;
                case(6)
                    obj.Plot = plot;
                    obj.Valmin = minValue;
                    obj.Valmax = maxValue;
                case(7)
                    obj.Plot = plot;
                    obj.Valmin = minValue;
                    obj.Valmax = maxValue;
                    obj.Value = init .* ones(row, col);
            end
        end
        function update(obj, delta)
            obj.Value = obj.Value + delta;
        end
    end
end

