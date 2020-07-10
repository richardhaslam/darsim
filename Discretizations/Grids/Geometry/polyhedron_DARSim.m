% Class of polyhedron for DARSim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DARSim Reservoir Simulator
% Author: Mousa HosseiniMehr
% TU Delft
% Created: 2020-03-09
% Last modified: 2020-06-17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef polyhedron_DARSim < handle
    properties
        NumOfVertex
        NumOfFace
    end
    methods
        %%
        function obj = polyhedron_DARSim(n_vertex,n_face)
            if nargin >= 2
                obj.InitializePolyhedron(n_vertex,n_face);
            end
        end
        %%
        function InitializePolyhedron(obj,n_vertex,n_face)
            obj.NumOfVertex = n_vertex;
            obj.NumOfFace   = n_face;
        end
        %%
        function [Geostatus, IntersectPoints] = Obtain_Polyhedron_LineSegment_Intersection(obj,LineSegment,Epsilon)
            % The algorithm used in this section is implemented from the method found in this article:
            % http://geomalgorithms.com/a13-_intersect-4.html
            
            % Initializing the output variables
            Geostatus.haveIntersect = NaN;
            IntersectPoints         = [];
            
            tE = 0; % for the maximum entering segment parameter
            tL = 1; % for the minimum leaving segment parameter
            for i= 1 : obj.NumOfFace
                % Checking if the normal vector is pointing outward, and if not, fixing it:
                if dot( obj.Centroid - obj.Face(i).PointM , obj.Face(i).nVec) > 0
                    obj.Face(i).nVec = - obj.Face(i).nVec;
                end
                
                N = - dot( LineSegment.PointA - obj.Face(i).PointA , obj.Face(i).nVec);
                D =   dot( LineSegment.AB_vec , obj.Face(i).nVec );
                
                if abs(D) < Epsilon
                    % Then, the lineSegment is parallel to this face
                    if N < 0
                        % Then, PointA of the lineSegment is outside of this face
                        % Getting out of the entire function for this hexahedron as
                        % this lineSegment will not intersect this hexahedron.
                        Geostatus.haveIntersect = 0;
                        return;
                    else
                        % Then, the llineSegment cannot enter or leave this
                        % face. Ignoring the rest of the loop for this face
                        % and moving to the next face.
                        continue;
                    end
                end
                
                t = N / D;
                if D < 0
                    % Then the lineSegment is entering the hexahedron across this face.
                    tE = max(tE, t);
                    if tE > tL
                        % Then, the lineSegment enters hexahedron after leaving.
                        % Getting out of the entire function for this hexahedron as
                        % this lineSegment will not intersect this hexahedron.
                        Geostatus.haveIntersect = 0;
                        return;
                    end
                else % D > 0
                    % Then the lineSegment is leaving the hexahedron across this face.
                    tL = min(tL, t);
                    if tL < tE
                        % Then, the lineSegment leaves hexahedron before entering.
                        % Getting out of the entire function for this hexahedron as
                        % this lineSegment will not intersect this hexahedron.
                        Geostatus.haveIntersect = 0;
                        return;
                    end
                end
            end
            
            if tE <= tL
                Geostatus.haveIntersect = 1;
                % We have intersection. The valid intersection points are:
                PointE = LineSegment.PointA + tE * LineSegment.AB_vec;
                PointL = LineSegment.PointA + tL * LineSegment.AB_vec;
                IntersectPoints = [PointE, PointL];
            end
        end
    end
end

