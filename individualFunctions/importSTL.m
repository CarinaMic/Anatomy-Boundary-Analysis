% Import stl files

% Input:    filename
%           Test before (f.e. in Blender) if the stl-file is correct and the normals show in the right direction (outwards).

% Output:   vertices
%           faces structure of the stl
%           normals of the faces
% Mirroring, if necessary

% Developed by C.Micheler,
% Department of Orthopaedics and Sportorthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich

function [vertices,faces,normals] = importSTL(filename,mirror)
    
    % Use the built-in MATLAB stlread function
    % Vertices/faces do not require restructuring (stlStruct function) with built-in stlread
    TR = stlread(filename);
    vertices = TR.Points;
    faces = TR.ConnectivityList;
    % Compute the face normals
    V1 = TR.Points(TR.ConnectivityList(:,1),:);
    V2 = TR.Points(TR.ConnectivityList(:,2),:);
    V3 = TR.Points(TR.ConnectivityList(:,3),:);
    Normals = cross(V2 - V1, V3 - V1, 2);
    % Normalize the normals
    normals = Normals ./ vecnorm(Normals,2,2);
    
    % Mirroring at yz-plane -> x-Coordinate * (-1)
    if mirror==1
        vertices(:,1) = (-1) * vertices(:,1);
        % Mesh: normal orientation to the center -> defines the outside of the mesh
        % The orientation of the normal is defined by the ordering of the vertices in a triangle (faces strcuture).
        % The ordering follows the right-hand rule.
        faces = TR.ConnectivityList(:,[1,3,2]);
        % Compute the face normals
        v1 = vertices(faces(:,1),:);
        v2 = vertices(faces(:,2),:);
        v3 = vertices(faces(:,3),:);
        normals = cross(v2 - v1, v3 - v1, 2);
        % Normalize the normals
        normals = normals ./ vecnorm(normals,2,2);
    end
    
    disp('STL data loaded')
end