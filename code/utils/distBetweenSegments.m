function min_dist = distBetweenSegments(P1, P2, P3, P4)
% DISTBETWEENSEGMENTS - Computes the minimum Euclidean distance between two 2D segments.
% 
% Syntax:
%   min_dist = distBetweenSegments(P1, P2, P3, P4)
%
% Inputs:
%   P1, P2 - Endpoints of the first segment [x, y] (e.g., Fault A Start/End)
%   P3, P4 - Endpoints of the second segment [x, y] (e.g., Fault B Start/End)
%
% Output:
%   min_dist - The shortest distance between the two segments. Returns 0 if they intersect.
%
% Logic:
%   1. Check for intersection using cross products.
%   2. If no intersection, the minimum distance must be between one of the endpoints
%      of one segment and the other segment.

    % --- 1. Intersection Check ---
    % If the segments cross each other, the distance is by definition zero.
    if segmentsIntersect(P1, P2, P3, P4)
        min_dist = 0;
        return;
    end

    % --- 2. Endpoint-to-Segment Distance ---
    % Calculate the distance from each of the 4 endpoints to the opposing segment.
    % The minimum of these 4 values is the overall minimum distance.
    dists = [
        pointToSegmentDist(P1, P3, P4), ... % Dist from P1 to segment P3-P4
        pointToSegmentDist(P2, P3, P4), ... % Dist from P2 to segment P3-P4
        pointToSegmentDist(P3, P1, P2), ... % Dist from P3 to segment P1-P2
        pointToSegmentDist(P4, P1, P2)      % Dist from P4 to segment P1-P2
    ];

    min_dist = min(dists);
end

% =========================================================
% HELPER FUNCTIONS
% =========================================================

function d = pointToSegmentDist(P, A, B)
% Calculates the shortest distance between point P and segment AB.
    v = B - A;
    w = P - A;
    
    % Projection of P onto the line containing AB
    c1 = dot(w, v);
    if c1 <= 0
        % Point is closest to endpoint A
        d = norm(P - A);
        return;
    end
    
    c2 = dot(v, v);
    if c2 <= c1
        % Point is closest to endpoint B
        d = norm(P - B);
        return;
    end
    
    % Point projection falls within the segment AB
    b = c1 / c2;
    Pb = A + b * v;
    d = norm(P - Pb);
end

function intersect = segmentsIntersect(P1, P2, P3, P4)
% Checks if segment P1-P2 and segment P3-P4 intersect.
% Uses the cross product orientation test.
    
    function cp = cross_product(a, b, c)
        % Standard 2D cross product of vectors (b-a) and (c-a)
        cp = (b(1)-a(1))*(c(2)-a(2)) - (b(2)-a(2))*(c(1)-a(1));
    end
    
    % Check orientations of points relative to segments
    cp1 = cross_product(P1, P2, P3);
    cp2 = cross_product(P1, P2, P4);
    cp3 = cross_product(P3, P4, P1);
    cp4 = cross_product(P3, P4, P2);
    
    % If points of one segment lie on opposite sides of the other segment, they intersect
    intersect = ((cp1*cp2 < 0) && (cp3*cp4 < 0));
end