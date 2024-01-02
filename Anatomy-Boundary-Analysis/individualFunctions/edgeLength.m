% Mean edge length of stl

% Input:    vertices
%           faces structure of the stl

% Output:   meanEdgeLength
%           stdEdgeLength

% Developed by C.Micheler,
% Department of Orthopaedics and Sportorthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich


function [meanEdgeLength,stdEdgeLength] = edgeLength(vertices,faces)
        % Calculate the edge length of the triangular faces in the mesh
        
        edgeLengths = []; % Initialize an array to store edge lengths
        for i = 1:size(faces,1)
            % Vertices of the triangle
            V1 = vertices(faces(i,1),:);
            V2 = vertices(faces(i,2),:);
            V3 = vertices(faces(i,3),:);
        
            % Calculate lengths of each edge of the triangle
            edgeLengths(end+1) = norm(V1 - V2); % Edge 1-2
            edgeLengths(end+1) = norm(V2 - V3); % Edge 2-3
            edgeLengths(end+1) = norm(V3 - V1); % Edge 3-1
        end
        
        meanEdgeLength = mean(edgeLengths);
        stdEdgeLength = std(edgeLengths);
        
        disp('Mean edge length calculated');

end