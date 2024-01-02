% Distances between the landmark points

% Input:    importData: struct with all landmarks

% Output:   vectors: Vectors of the landmarks to each other
%           distances: Distances between the landmark points
%           combis: Possible combinations of all landmarks

% Developed by C.Micheler,
% Department of Orthopaedics and Sportorthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich


function [vectors, distances, combis] = distanceLM(importData)
    % Define all possible landmarks
    allLandmarks = {'asis', 'tb', 'psis', 'piis', 'si', 'acentre', 'aiis', 'iim', ...
        'iimin', 'spp', 'spd', 'ti', 'fo', 'ci'};
    
    % Initialize landmark data with NaN
    landmarkData = NaN(length(allLandmarks), 3);
    
    % Populate available landmark data
    for idx = 1:length(allLandmarks)
        if isfield(importData, allLandmarks{idx})
            landmarkData(idx, :) = importData.(allLandmarks{idx});
        end
    end
    
    % Landmark combinations
    LMcount = length(allLandmarks);
    num = nchoosek(LMcount, 2); % Combination count based on available landmarks
    combis = nchoosek(1:LMcount, 2); % Combinations
    
    vectors = zeros(num, 3);
    distances = zeros(num, 1);
    
    for j = 1:num
        % Ensure both landmarks are available before calculating distance
        if ~any(isnan(landmarkData(combis(j, :), :)))
            % Landmark vectors
            vectors(j, :) = landmarkData(combis(j, 1), :) - ...
                landmarkData(combis(j, 2), :);
    
            % Euclidean distance (norm)
            distances(j) = norm(vectors(j, :));
        else
            vectors(j, :) = [NaN, NaN, NaN];
            distances(j) = NaN;
        end
    end
    
    disp('Landmark distances calculated');
end