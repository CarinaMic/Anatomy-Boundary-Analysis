% Distances between the landmarks and the centre of the acetabulum

% Input:    importData: struct with all landmarks

% Output:   vectors: Vectors of the landmarks to each other
%           distances: Distances between the landmark points
%           closestLM: Five closest landmarks to 'acentre'

% Developed by C.Micheler,
% Department of Orthopaedics and Sportorthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich

% Distances landmarks to acentre
function [vectors,distances,closestLM] = distanceCentreLM(importData)
    
    % Centre landmark
    if isfield(importData, 'acentre')
        landmarkCentre(1, :) = importData.acentre;
    else
        error('Landmark "acentre" is not available in the processed data.');
    end
    
    % Landmarks excluding the center
    allLandmarks = {'asis', 'tb', 'psis', 'piis', 'si', 'aiis', 'iim', ...
        'iimin', 'spp', 'spd', 'ti', 'fo', 'ci'};
    
    % Initialize the matrix to store landmark data with NaN
    LMcount = length(allLandmarks); % Landmark combinations
    landmarkData = NaN(LMcount, 3);
    
    % Fetch available landmark data using a loop
    for idx = 1:LMcount
        if isfield(importData, allLandmarks{idx})
            landmarkData(idx, :) = importData.(allLandmarks{idx});
        end
    end
    
    % Landmark vectors
    vectors = landmarkCentre - landmarkData;
    
    % If any component of a landmark is NaN
    vectors(any(isnan(vectors), 2), :) = NaN;
    
    % Euclidean distance (norm)
    distances = sqrt(sum(vectors.^2, 2));
    
    % Finding the five closest landmarks to 'acentre'
    [sortedDistances, sortedIndices] = sort(distances);
    closestLandmarks = allLandmarks(sortedIndices(1:5)); % Top 5 closest landmarks
    
    % Store closest landmarks
    closestLM = closestLandmarks;
    
    disp(['Landmark distances to centre calculated: pelvis ', num2str(pelvisNum)]);

end