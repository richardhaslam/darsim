% Class of polyhedron for DARSim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim Reservoir Simulator
%Author: Mousa HosseiniMehr
%TU Delft
%Created: 2020-03-09
%Last modified: 2020-03-09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef polyhedron_DARSim < handle
    properties
        N_Vertex
        N_Face
    end
    methods
        function obj = polyhedron_DARSim(n_vertex,n_face)
            if nargin==2
                obj.InitializePolyhedron(n_vertex,n_face);
            end
        end
        function InitializePolyhedron(obj,n_vertex,n_face)
            obj.N_Vertex = n_vertex;
            obj.N_Face = n_face;
        end
    end
end

