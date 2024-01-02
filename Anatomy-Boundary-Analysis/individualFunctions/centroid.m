% Centroid of the faces (for the normal vector)

% Input:    vertices
%           faces structure of the stl

% Output:   centreFaces (centre of the triangle)

% Developed by C.Micheler,
% Department of Orthopaedics and Sportorthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich


function centreFaces = centroid(vertices,faces)
    % Calculate the centroid of the triangular faces in the mesh
    V1 = vertices(faces(:,1),:);
    V2 = vertices(faces(:,2),:);
    V3 = vertices(faces(:,3),:);
    faceVertices = cat(3, V1, V2, V3);
    centreFaces = mean(faceVertices, 3);
    
    disp('Centroids of the faces calculated');

end