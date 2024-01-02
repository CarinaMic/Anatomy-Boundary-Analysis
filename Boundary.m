classdef Boundary
    % Boundary class to build a boundary box/sphere around model or to calculate distances of landmarks
    % For scaling of the model
    
    % Developed by C.Micheler,
    % Department of Orthopaedics and Sportorthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
    % Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich
    
    
    properties 
        cuboid = struct();      % Boundary cuboid
        sphere = struct();      % Boundary sphere
        landmarks = struct();   % Landmarks distance  
    end
    
    methods
        %% Constructer: generate object
        function obj = Boundary()
            disp('class boundary initialized')
        end  
        
        %% Boundary box / cuboid: 
        function obj = boundBoxHull(obj,pelvisNum,points,level)
            
            % Boundary box
            % https://de.mathworks.com/matlabcentral/fileexchange/18264-minimal-bounding-box
            % https://blogs.mathworks.com/pick/2019/10/25/minimal-bounding-box/
            
            % metric - (OPTIONAL) - single letter character flag which
            %        denotes the use of minimal volume, surface or sum of edges as the
            %        metric to be minimized. metric may be either 'v', 's' or 'e',
            %        capitalization is ignored.
            %        DEFAULT: 'v'    ('volume')
            % level - (OPTIONAL) - either 1, 2, 3 or 4. This denotes the level of
            %       reliability of the resulting minimal bounding box.
            %       '1' denotes the search for a bounding box, where one side of the
            %           box coincides with a face of the convex hull. (fast)
            %       '2' like '1', but also incorporates a search through each pair of edges
            %           which form a plane which coincides with one side of the box (slow)
            %       '3' like '1', but also incorporates a search through each edge which is
            %           parallel to an edge of the box (reasonably fast)
            %       '4' like '1', '2' and '3' together. (slowest) (Never needed that.)
            %       It depends on the application, what should be chosen here.
            %       See the example at the end for the effects of this parameter.
            %       DEFAULT: '3'
            [obj.cuboid.hull.rotmat,obj.cuboid.hull.cornerpoints,obj.cuboid.hull.vol,...
                obj.cuboid.hull.surf,obj.cuboid.hull.sumAllEdges] = ...
                minboundbox(points(:,1),points(:,2),points(:,3),'v',level);
            
            % Cuboid parameters
            [obj.cuboid.hull.volume, obj.cuboid.hull.surface, ...
            obj.cuboid.hull.edgeVector, obj.cuboid.hull.edgeLength,...
                obj.cuboid.hull.sumEdges, obj.cuboid.hull.diagonal] = paraCuboid(obj.cuboid.hull.cornerpoints); 
            
            % Triangulation of the cornerpoints (convhull: matlab function)
            [obj.cuboid.hull.tri,obj.cuboid.hull.vol] = convhull(obj.cuboid.hull.cornerpoints(:,1),...
                obj.cuboid.hull.cornerpoints(:,2), obj.cuboid.hull.cornerpoints(:,3),'Simplify',true);
            
            disp(['boundary box (hull) calculated: pelvis ',num2str(pelvisNum)])
            
        end
        
        %% Boundary box / cuboid 2: singular value decomposition (SVD) of the points on the convex hull 
        function obj = boundBoxSVD(obj,pelvisNum,points)
        
            % https://de.mathworks.com/matlabcentral/fileexchange/64417-calc_oriboundingbox-data?s_tid=srchtitle
            obj.cuboid.SVD.cornerpoints = calc_OriBoundingBox(points);
            
            % Cuboid parameters
            [obj.cuboid.SVD.volume, obj.cuboid.SVD.surface,...
                obj.cuboid.SVD.edgeVector, obj.cuboid.SVD.edgeLength,...
                obj.cuboid.SVD.sumEdges, obj.cuboid.SVD.diagonal] = paraCuboid(obj.cuboid.SVD.cornerpoints); 
            
            % Triangulation of the cornerpoints (convhull: matlab function)
            [obj.cuboid.SVD.tri,obj.cuboid.SVD.vol] = convhull(obj.cuboid.SVD.cornerpoints(:,1),...
                obj.cuboid.SVD.cornerpoints(:,2), obj.cuboid.SVD.cornerpoints(:,3),'Simplify',true);
                        
            disp(['boundary box (SVD) calculated: pelvis ',num2str(pelvisNum)])
            
        end
               
        %% Boundary sphere 
        function obj = boundSphere(obj,pelvisNum,coordinates)
         
            % Minimal boundary sphere
            % https://de.mathworks.com/matlabcentral/fileexchange/48725-exact-minimum-bounding-spheres-and-circles
            [obj.sphere.R,obj.sphere.centre,~] = ExactMinBoundSphere3D(coordinates);
            
            disp(['boundary sphere calculated: pelvis ',num2str(pelvisNum)])
            
        end
        
         %% Distances landmarks
        function obj = distanceLM(obj, pelvisNum, processedData)
            % Define all possible landmarks
            allLandmarks = {'asis', 'tb', 'psis', 'piis', 'si', 'acentre', 'aiis', 'iim', ...
                'iimin', 'spp', 'spd', 'ti', 'fo', 'ci'};

            % Initialize landmark data with NaN
            landmarkData = NaN(length(allLandmarks), 3);

            % Populate available landmark data
            for idx = 1:length(allLandmarks)
                if isfield(processedData, allLandmarks{idx})
                    landmarkData(idx, :) = processedData.(allLandmarks{idx});
                end
            end

            % Landmark combinations
            LMcount = length(allLandmarks);
            num = nchoosek(LMcount, 2); % Combination count based on available landmarks
            obj.landmarks.all.combis = nchoosek(1:LMcount, 2); % Combinations

            vectors = zeros(num, 3);
            distances = zeros(num, 1);

            parfor j = 1:num
                % Ensure both landmarks are available before calculating distance
                if ~any(isnan(landmarkData(obj.landmarks.all.combis(j, :), :)))
                    % Landmark vectors
                    vectors(j, :) = landmarkData(obj.landmarks.all.combis(j, 1), :) - ...
                        landmarkData(obj.landmarks.all.combis(j, 2), :);

                    % Euclidean distance (norm)
                    distances(j) = norm(vectors(j, :));
                else
                    vectors(j, :) = [NaN, NaN, NaN];
                    distances(j) = NaN;
                end
            end

            obj.landmarks.all.distance = distances;
            obj.landmarks.all.vectors = vectors;

            disp(['Landmark distances calculated: pelvis ', num2str(pelvisNum)]);
        end
        
        %% Distances landmarks to acentre
        function obj = distanceCentreLM(obj,pelvisNum,processedData)

            % Centre landmark
            if isfield(processedData, 'acentre')
                landmarkCentre(1, :) = processedData.acentre;
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
                if isfield(processedData, allLandmarks{idx})
                    landmarkData(idx, :) = processedData.(allLandmarks{idx});
                end
            end

            % Landmark vectors
            vectors = landmarkCentre - landmarkData;

            % If any component of a landmark is NaN
            vectors(any(isnan(vectors), 2), :) = NaN;

            % Euclidean distance (norm)
            distances = sqrt(sum(vectors.^2, 2));

            obj.landmarks.centre.distance = distances;
            obj.landmarks.centre.vectors = vectors;

            % Finding the five closest landmarks to 'acentre'
            [sortedDistances, sortedIndices] = sort(distances);
            closestLandmarks = allLandmarks(sortedIndices(1:5)); % Top 5 closest landmarks

            % Store closest landmarks
            obj.landmarks.centre.closestLM = closestLandmarks;

            disp(['Landmark distances to centre calculated: pelvis ', num2str(pelvisNum)]);
            
        end
        
    end
    
end