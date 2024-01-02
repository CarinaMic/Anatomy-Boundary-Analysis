% Calculating the parameters of the cuboid

% Developed by C.Micheler,
% Department of Orthopaedics and Sportorthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich
    
% Cuboid parameters
function [volume, surface, edgeVector, edgeLength, sumEdges, diagonal] = paraCuboid(cornerpoints)

    % % Unsorted cornerpoints
    % % Vectors in cuboid; (1,:) = start point
    % for j = 2:8 % 7 further cornerpoints
    %     % Cornerpoints vectors
    %     vectorsLower(j-1,:) = cornerpoints(j,:) - cornerpoints(1,:);
    %     % Euclidean distance (norm)
    %     distances(j-1,:) = norm(vectorsLower(j-1,:));
    % end
    % % Minimal vector length -> first edge vector
    % [minDistance, indexMinDistance] = min(distances);
    % edgeVector(1,:) = vectorsLower(indexMinDistance,:);
    % edgeLength(1,:) = distances(indexMinDistance,:);
    % % Two vectors are perpendicular to each other if their scalar product is zero
    % for j = 1:7
    %     scalar(j,:) = dot(edgeVector(1,:),vectorsLower(j,:));
    %     CosTheta(j,:) = max(min(dot(edgeVector(1,:),vectorsLower(j,:))/...
    %         (norm(edgeVector(1,:))*norm(vectorsLower(j,:))),1),-1);
    %     ThetaInDegrees(j,:) = real(acosd(CosTheta(j,:))); % angle of the vectors to each other
    % end
    % 
    % % Three vectors perpendicular to first edge vector
    % [rowEdges, ~] = find(ThetaInDegrees>89 & ThetaInDegrees<91);
    % edge(1,:) = distances(rowEdges(1),:);
    % edge(2,:) = distances(rowEdges(2),:);
    % edge(3,:) = distances(rowEdges(3),:);
    % % The vector with the greatest distance is the diagonal of the other two vectors
    % [maxEdge,indexMaxEdge] = max(edge);
    % rowEdges(indexMaxEdge) = [];
    % 
    % % Second and third edge vector of length, width and height
    % edgeVector(2,:) = vectorsLower(rowEdges(1),:);
    % edgeVector(3,:) = vectorsLower(rowEdges(2),:);
    % edgeLength(2,:) = distances(rowEdges(1),:);
    % edgeLength(3,:) = distances(rowEdges(2),:);

    
    % Three sorted edges (sorted cornerpoints)
    % Sorted cornerpoints in local COSY (box):
    % xmin,ymin,zmin
    % xmax,ymin,zmin
    % xmax,ymax,zmin
    % xmin,ymax,zmin
    % xmin,ymin,zmax
    % xmax,ymin,zmax
    % xmax,ymax,zmax
    % xmin,ymax,zmax
    edgeVector(1,:) = cornerpoints(2,:) - cornerpoints(1,:); % x
    edgeVector(2,:) = cornerpoints(4,:) - cornerpoints(1,:); % y
    edgeVector(3,:) = cornerpoints(5,:) - cornerpoints(1,:); % z
    edgeLength(1,:) = norm(edgeVector(1,:));
    edgeLength(2,:) = norm(edgeVector(2,:));
    edgeLength(3,:) = norm(edgeVector(3,:));

    % Sum of edges
    sumEdges = edgeLength(1,:) + edgeLength(2,:) + edgeLength(3,:);
    % Surface
    surface = 2*edgeLength(1,:)*edgeLength(2,:) + 2*edgeLength(1,:)*edgeLength(3,:) + ...
        2*edgeLength(2,:)*edgeLength(3,:);
    % Volume
    volume = edgeLength(1,:) * edgeLength(2,:) * edgeLength(3,:);
    % Diagonal: triangle of first edge Vector (=min) and vector with the greatest distance (maxEdge)
    % diagonal = sqrt(edgeLength(1,:)^2 + maxEdge^2);
    diagonal = norm(cornerpoints(7,:) - cornerpoints(1,:)); % sorted cornerpoints)
    
end