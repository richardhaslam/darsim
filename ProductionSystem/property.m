classdef property < handle
    % This is a generic Property
    properties
        Value
    end
    methods
        function obj = property(row, col)
            % Creates matrix with property values
            obj.Value = zeros(row, col);
        end
    end
end

