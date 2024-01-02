%% Matlab initialisation script

% Project to analyse the bounding methods for anatomical structures (e.g. here on the pelvis)

% Developed by C.Micheler,
% Department of Orthopaedics and Sportorthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich


% Clear environment and close figures
clc;
clear variables;
close all;

% Add necessary paths
addpath('./Functions/'); % Add Function folder

% Parallel computing setup
if isempty(gcp('nocreate'))
    parpool('local'); % Start parallel pool if not already running
end

% Define TUM color palette
TUMcolors = struct( ...
    'blue300', [0, 101/255, 189/255], ...
    'blue540', [0, 51/255, 89/255], ...
    'blue301', [0, 82/255, 147/255], ...
    'blue542', [100/255, 160/255, 200/255], ...
    'blue283', [152/255, 198/255, 234/255], ...
    'grey80',  [51/255, 51/255, 51/255], ...
    'grey50',  [128/255, 128/255, 128/255], ...
    'grey20',  [204/255, 204/255, 204/255], ...
    'cream',   [218/255, 215/255, 203/255], ...
    'orange',  [227/255, 114/255, 34/255], ...
    'green',   [162/255, 173/255, 0/255] ...
);

%% %%%%%%%%%%%% Class Pelvis %%%%%%%%%%%%%%%

% Folder of the geometries
filedir = './GeometriesPelves/';
% Lists of stl and xlsx files
stlfiles = dir(fullfile(filedir, '*.stl'));
xlsfiles = dir(fullfile(filedir, '*.xlsx'));
% Use stl/xls files count to determine data count
dataCount = length(stlfiles);
%dataCount = length(xlsfiles);

% Check for mismatch or absence of files
if isempty(stlfiles)
    error('No stl files found in the specified directory.');
elseif isempty(xlsfiles)
    error('No xlsx files found in the specified directory.');
elseif length(stlfiles) ~= length(xlsfiles)
    error('Number of stl files does not match number of xlsx files.');
end

% Initalization: Generate object array
pelvis(dataCount,1) = Pelvis;

%% %%%%%%%%%%%% Class Import %%%%%%%%%%%%%%%
% Import pelvis data

% File paths
filePath_pelvisData = fullfile(filedir, {xlsfiles.name}); % Pelvis data/infos
filePath_pelvis = fullfile(filedir, {stlfiles.name}); % Pelvis geometry data

% Parallel loop for data import and processing
parfor i = 1:dataCount

   % Initalization: Generate object
   pelvisImport = Import(i); % cache (to ensure that the parallel calculations are independent) 
      
   % Import the pelvis data/infos
   pelvisImport = pelvisImport.importData(i,filePath_pelvisData{i});
  
   % Import the pelvis geometry data
   pelvisImport = pelvisImport.importSTL(i, filePath_pelvis{i}) ...
                               .centroid(i) ... % Centroid of the faces (for the normal vector)
                               .edgeLength(i); % Edge length of  mesh

   % Mean edge length
   allEdgeLength(i) = pelvisImport.processed.meanEdgeLength;

   % Store the result back in pelvis
   pelvis(i).import = pelvisImport;

end

% Mean edge Length
allMeanEdge = mean(allEdgeLength);
allStdEdge = std(allEdgeLength);

%% Display the pelvis with the face normals (optional), acetabulum centre and pelvis landmarks (for control)

% List of landmarks with their respective colors
landmarkList = {'acentre', 'asis', 'tb', 'psis', 'piis', 'si', ...
                'aiis', 'iim', 'iimin', 'spp', 'spd', 'ti', ...
                'fo', 'ci'};
landmarkColors = {TUMcolors.blue300, TUMcolors.orange, TUMcolors.orange, ...
                  TUMcolors.orange, TUMcolors.orange, TUMcolors.orange, ...
                  TUMcolors.orange, TUMcolors.orange, TUMcolors.orange, ...
                  TUMcolors.orange, TUMcolors.orange, TUMcolors.orange, ...
                  TUMcolors.orange, TUMcolors.orange};

% List of auxiliary landmarks with yellow color
auxLandmarkList = {'fo_p', 'fo_i', 'fo_a', 'fo_ip'};
auxLandmarkColors = {TUMcolors.green, TUMcolors.green, ...
                     TUMcolors.green, TUMcolors.green};

parfor i = 1:dataCount
    figure('Visible','off')
   
    % Pelvis
    patch('Faces', pelvis(i).import.processed.faces, ...
          'Vertices', pelvis(i).import.processed.vertices, ...
          'FaceColor', [0.9 0.75 0.68], ...
          'FaceAlpha', 1, ...
          'EdgeColor', TUMcolors.grey50, ...
          'EdgeAlpha', 0.25);
    
    hold on;
    
    % Plotting main landmarks
    for j = 1:length(landmarkList)
        landmark = pelvis(i).import.processed.(landmarkList{j});
        plot3(landmark(1), landmark(2), landmark(3), '*', 'Color', landmarkColors{j});
    end
    
    % Plotting auxiliary landmarks
    for j = 1:length(auxLandmarkList)
        auxLandmark = pelvis(i).import.processed.(auxLandmarkList{j});
        plot3(auxLandmark(1), auxLandmark(2), auxLandmark(3), '*', 'Color', auxLandmarkColors{j});
    end
    
    % Format and display properties
    xlabel('x')
    ylabel('y')
    zlabel('z')
    grid on
    daspect([1, 1, 1]); % Equal aspect ratio for the axes
    view(3);
    hold off;
    
    % Save figure (figure invisible, but saved visible)
    set(gcf, 'Visible', 'off', 'CreateFcn', 'set(gcf,''Visible'',''on'')') 
    savefig(['./Figures/Pelvis(', num2str(i), ').fig'])
end

clear landmarkList landmarkColors auxLandmarkList auxLandmarkColors landmark auxLandmark

%% %%%%%%%%%%% Class Boundary %%%%%%%%%%%%%%
%% Boundary cuboid: based on convex hull and its faces or edges

% Bounding box around pelves
% level - (last input) - either 1, 2, 3 or 4. This denotes the level of
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
level = 1;

parfor i = 1:dataCount

    pelvis(i).boundaries = pelvis(i).boundaries.boundBoxHull(i, pelvis(i).import.processed.vertices,level); 
   
end

%% Display the pelvis with the boundary box: based on convex hull and its faces or edges (for control)
% Box / cuboid
parfor i = 1 : dataCount
    figure('Visible','off')   
    patch('Faces', pelvis(i).import.processed.faces, ...
          'Vertices', pelvis(i).import.processed.vertices, ...
          'FaceColor', [0.9 0.75 0.68], ...
          'FaceAlpha', 1, ...
          'EdgeColor', TUMcolors.grey50, ...
          'EdgeAlpha', 0.25);
    hold on

    % Display bounding box
    trisurf(pelvis(i).boundaries.cuboid.hull.tri,...
        pelvis(i).boundaries.cuboid.hull.cornerpoints(:,1),...
        pelvis(i).boundaries.cuboid.hull.cornerpoints(:,2),...
        pelvis(i).boundaries.cuboid.hull.cornerpoints(:,3),...
        'FaceColor',TUMcolors.blue300,'EdgeColor',TUMcolors.blue300,'FaceAlpha',0.25);

    % Display cuboid edges
    for j=1:3
        quiver3(pelvis(i).boundaries.cuboid.hull.cornerpoints(1,1),...
            pelvis(i).boundaries.cuboid.hull.cornerpoints(1,2),...
            pelvis(i).boundaries.cuboid.hull.cornerpoints(1,3),... % start point
            pelvis(i).boundaries.cuboid.hull.edgeVector(j,1),...
            pelvis(i).boundaries.cuboid.hull.edgeVector(j,2),...
            pelvis(i).boundaries.cuboid.hull.edgeVector(j,3),0,'LineWidth',2)
    end

    % Format and display properties
    xlabel('x')
    ylabel('y')
    zlabel('z')
    grid on
    daspect([1, 1, 1]); % Equal aspect ratio for the axes
    view(3);
    hold off;

    % Save figure (figure unvisible, but saved visible)
    set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')') 
    savefig(['./Figures/Pelvis(',num2str(i),')CuboidMinHull.fig'])
end

%% Picture BMT paper 2022
% i = 28; %X00110
% figure
% patch('Faces',pelvis(i).import.processed.faces,...
%     'Vertices',pelvis(i).import.processed.vertices,...
%     'FaceColor',[0.8 0.8 0.8], ...      % Face color
%     'FaceAlpha',1,...                   % Transparency of the faces
%     'EdgeColor',[0.6 0.6 0.6],...       % Edge color
%     'EdgeAlpha',0.5);                   % Transparency of the edges
% hold on
% 
% % Display bounding box
% trisurf(pelvis(i).boundaries.cuboid.hull.tri,...
%     pelvis(i).boundaries.cuboid.hull.cornerpoints(:,1),...
%     pelvis(i).boundaries.cuboid.hull.cornerpoints(:,2),...
%     pelvis(i).boundaries.cuboid.hull.cornerpoints(:,3),...
%     'FaceColor',[0 0.396 0.729],'EdgeColor',[0 0.396 0.729],'FaceAlpha',0.2);
% 
% % Format and display properties
% xlabel('x')
% ylabel('y')
% zlabel('z')
% grid off
% daspect([1, 1, 1]); % Equal aspect ratio for the axes
% view(3);
% hold off;

%% Boundary cuboid: singular value decomposition (SVD) of the points on the convex hull
parfor i = 1:dataCount
    
    % Bounding box around pelves
    pelvis(i).boundaries = pelvis(i).boundaries.boundBoxSVD(i,pelvis(i).import.processed.vertices);
    
end

%% Display the pelvis with the boundary box: SVD (for control) 
% Box / cuboid
parfor i = 1:dataCount
    figure('Visible','off')
    patch('Faces',pelvis(i).import.processed.faces,...
        'Vertices',pelvis(i).import.processed.vertices,...
        'FaceColor', [0.9 0.75 0.68], ...
        'FaceAlpha', 1, ...
        'EdgeColor', TUMcolors.grey50, ...
        'EdgeAlpha', 0.25);
    hold on

    % Display bounding box
    trisurf(pelvis(i).boundaries.cuboid.SVD.tri,...
        pelvis(i).boundaries.cuboid.SVD.cornerpoints(:,1),...
        pelvis(i).boundaries.cuboid.SVD.cornerpoints(:,2),...
        pelvis(i).boundaries.cuboid.SVD.cornerpoints(:,3),...
        'FaceColor',TUMcolors.blue300,'EdgeColor',TUMcolors.blue300,'FaceAlpha',0.25);

    % Display cuboid edges
    for j=1:3
        quiver3(pelvis(i).boundaries.cuboid.SVD.cornerpoints(1,1),...
            pelvis(i).boundaries.cuboid.SVD.cornerpoints(1,2),...
            pelvis(i).boundaries.cuboid.SVD.cornerpoints(1,3),... % start point
            pelvis(i).boundaries.cuboid.SVD.edgeVector(j,1),...
            pelvis(i).boundaries.cuboid.SVD.edgeVector(j,2),...
            pelvis(i).boundaries.cuboid.SVD.edgeVector(j,3),0,'LineWidth',2)
    end
    
    % Format and display properties
    xlabel('x')
    ylabel('y')
    zlabel('z')
    grid on
    daspect([1, 1, 1]); % Equal aspect ratio for the axes
    view(3);
    hold off;

    % Save figure (figure unvisible, but saved visible)
    set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')') 
    savefig(['./Figures/Pelvis(',num2str(i),')CuboidSVD.fig'])
end

%% Picture BMT paper 2022
% i = 13; %180134
% figure
% patch('Faces',pelvis(i).import.processed.faces,...
%     'Vertices',pelvis(i).import.processed.vertices,...
%     'FaceColor',[0.8 0.8 0.8], ...      % Face color
%     'FaceAlpha',1,...                   % Transparency of the faces
%     'EdgeColor',[0.6 0.6 0.6],...       % Edge color
%     'EdgeAlpha',0.5);                   % Transparency of the edges
% hold on
% 
% % Display bounding box
% trisurf(pelvis(i).boundaries.cuboid.SVD.tri,...
%     pelvis(i).boundaries.cuboid.SVD.cornerpoints(:,1),...
%     pelvis(i).boundaries.cuboid.SVD.cornerpoints(:,2),...
%     pelvis(i).boundaries.cuboid.SVD.cornerpoints(:,3),...
%     'FaceColor',[0 0.396 0.729],'EdgeColor',[0 0.396 0.729],'FaceAlpha',0.2);
% 
% % Format and display properties
% xlabel('x')
% ylabel('y')
% zlabel('z')
% grid off
% daspect([1, 1, 1]); % Equal aspect ratio for the axes
% view(3);
% hold off;

%% Boundary sphere
parfor i = 1:dataCount
   
    % Minimal bounding sphere around pelves
    pelvis(i).boundaries = pelvis(i).boundaries.boundSphere(i,pelvis(i).import.processed.vertices);
     
end

%% Display the pelvis with the boundary sphere (for control)
% Sphere

parfor i = 1:dataCount
    figure('Visible','off')
    patch('Faces',pelvis(i).import.processed.faces,...
        'Vertices',pelvis(i).import.processed.vertices,...
        'FaceColor',[0.8 0.8 0.8], ...      % Face color
        'FaceAlpha',1,...                   % Transparency of the faces
        'EdgeColor',[0.6 0.6 0.6],...       % Edge color
        'EdgeAlpha',0.5);                   % Transparency of the edges
    hold on
    [x,y,z] = sphere;
    % Scale to desire radius
    x = x * pelvis(i).boundaries.sphere.R;
    y = y * pelvis(i).boundaries.sphere.R;
    z = z * pelvis(i).boundaries.sphere.R;
    % C: Translate sphere to new location
    % Figure: Plot as surface
    surf(x+pelvis(i).boundaries.sphere.centre(1),y+pelvis(i).boundaries.sphere.centre(2),...
        z+pelvis(i).boundaries.sphere.centre(3), ...
        'FaceColor',TUMcolors.blue300,'EdgeColor',TUMcolors.blue300,'FaceAlpha',0.25);
    hold on
    % scatter3(pelvis(i).import.processed.vertices(:,1),...
    %     pelvis(i).import.processed.vertices(:,2),...
    %     pelvis(i).import.processed.vertices(:,3));
    % hold on
    
    % Format and display properties
    xlabel('x')
    ylabel('y')
    zlabel('z')
    grid on
    axis equal;
    daspect([1, 1, 1]); % Equal aspect ratio for the axes
    view(3);
    hold off;

    % Save figure (figure unvisible, but saved visible)
    set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')') 
    savefig(['./Figures/Pelvis(',num2str(i),')Sphere.fig'])
end

%% Picture BMT paper 2022
% i=34; %Y00199
% figure
% patch('Faces',pelvis(i).import.processed.faces,...
%     'Vertices',pelvis(i).import.processed.vertices,...
%     'FaceColor',[0.8 0.8 0.8], ...      % Face color
%     'FaceAlpha',1,...                   % Transparency of the faces
%     'EdgeColor',[0.6 0.6 0.6],...       % Edge color
%     'EdgeAlpha',0.5);                   % Transparency of the edges
% hold on
% [x,y,z] = sphere;
% % Scale to desire radius
% x = x * pelvis(i).boundaries.sphere.R;
% y = y * pelvis(i).boundaries.sphere.R;
% z = z * pelvis(i).boundaries.sphere.R;
% % C: Translate sphere to new location
% % Figure: Plot as surface
% surf(x+pelvis(i).boundaries.sphere.centre(1),y+pelvis(i).boundaries.sphere.centre(2),...
%     z+pelvis(i).boundaries.sphere.centre(3),...
%     'FaceColor',[0 0.396 0.729],'FaceAlpha', 0.2, 'EdgeColor',[0 0.396 0.729])
% 
% % Format and display properties
% xlabel('x')
% ylabel('y')
% zlabel('z')
% grid off
% axis equal;
% daspect([1, 1, 1]); % Equal aspect ratio for the axes
% view(3);
% hold off;

%% Distance landmarks
parfor i = 1:dataCount
   
    pelvis(i).boundaries = pelvis(i).boundaries.distanceLM(i,pelvis(i).import.processed);
    LMallDistance(:,i) = pelvis(i).boundaries.landmarks.all.distance(:,1)

end

%% Display the landmark distances in bar graph
parfor i = 1:dataCount
    figure('Visible','off')
    bar(1:91,pelvis(i).boundaries.landmarks.all.distance(:,1))

    % Format and display properties
    xlabel('Distance Nr.')
    ylabel('Distance in mm')
    grid on

    % Save figure (figure unvisible, but saved visible)
    set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')') 
    savefig(['./Figures/Pelvis(',num2str(i),')LM.fig'])
end

figure('Visible','off')
bar(1:91,LMallDistance)

% Format and display properties
xlabel('Distance Nr.')
ylabel('Distance in mm')
grid on

% Save figure (figure unvisible, but saved visible)
set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')')
savefig('./Figures/PelvisLM.fig')

%% Distance landmarks to acetabulum centre (acentre)
parfor i = 1:dataCount
   
    pelvis(i).boundaries = pelvis(i).boundaries.distanceCentreLM(i,pelvis(i).import.processed);
    LMcentreDistance(:,i) = pelvis(i).boundaries.landmarks.centre.distance(:,1)
   
end

% Counts for closest landmarks to acentre
allLandmarks = {'asis', 'tb', 'psis', 'piis', 'si', 'acentre', 'aiis', 'iim', ...
    'iimin', 'spp', 'spd', 'ti', 'fo', 'ci'};
for j = 1:length(allLandmarks)
    LMcounts.(allLandmarks{j}) = 0;
end
% methodBound: 'distanceCentreLM'
for i = 1:dataCount
    % Counts for closest landmarks to acentre
    closestLM = pelvis(i).boundaries.landmarks.centre.closestLM;
    for j = 1:length(closestLM)
        landmark = closestLM{j};
        if isfield(LMcounts, landmark)
            LMcounts.(landmark) = LMcounts.(landmark) + 1;
        end
    end
end

%% Display the landmark distance to acetabulum centre in bar graph
parfor i = 1:dataCount
    figure('Visible','off')
    bar(1:13,pelvis(i).boundaries.landmarks.centre.distance(:,1))

    % Format and display properties
    xlabel('Distance Nr.')
    ylabel('Distance in mm')
    grid on

    % Save figure (figure unvisible, but saved visible)
    set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')') 
    savefig(['./Figures/Pelvis(',num2str(i),')LMcentre.fig'])
end

figure('Visible','off')
bar(1:13,LMcentreDistance)

% Format and display properties
xlabel('Distance Nr.')
ylabel('Distance in mm')
grid on

% Save figure (figure unvisible, but saved visible)
set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')')
savefig('./Figures/PelvisLMcentre.fig')

%% Evaluation: Scaling factor / Size comparison to reference

% Reference
% Reference diagonal (bounding box): first pelvis data
refCuboidHullDiagonal = pelvis(1).boundaries.cuboid.hull.diagonal;
refCuboidSVDDiagonal = pelvis(1).boundaries.cuboid.SVD.diagonal;
% Reference radius (sphere): first pelvis data
refSphereR = pelvis(1).boundaries.sphere.R;
% Reference distances (landmarks): first pelvis data
refLMdistance = pelvis(1).boundaries.landmarks.all.distance;
refLMdistanceCentre = pelvis(1).boundaries.landmarks.centre.distance;

% Size comparison to reference
for i = 1:dataCount 
    
    % Cuboid/Bounding Box Hull
    pelvis(i).boundaries.cuboid.hull.size = pelvis(i).boundaries.cuboid.hull.diagonal / ...
        refCuboidHullDiagonal;                                                          % Size to reference (percent)
    pelvis(i).boundaries.cuboid.hull.scale = 1 / pelvis(i).boundaries.cuboid.hull.size; % Scaling factor
    pelvisScale.cuboidHull.diagonal(i,1) = pelvis(i).boundaries.cuboid.hull.diagonal;   % Cuboid diagonal
    pelvisScale.cuboidHull.size(i,1) = pelvis(i).boundaries.cuboid.hull.size;           % Size to reference
    pelvisScale.cuboidHull.scale(i,1) = pelvis(i).boundaries.cuboid.hull.scale;         % Scaling factor 
    
    % Cuboid/Bounding Box SVD
    pelvis(i).boundaries.cuboid.SVD.size = pelvis(i).boundaries.cuboid.SVD.diagonal / ...
        refCuboidSVDDiagonal;                                                           % Size to reference (percent)
    pelvis(i).boundaries.cuboid.SVD.scale = 1 / pelvis(i).boundaries.cuboid.SVD.size;   % Scaling factor
    pelvisScale.cuboidSVD.diagonal(i,1) = pelvis(i).boundaries.cuboid.SVD.diagonal;     % Cuboid diagonal
    pelvisScale.cuboidSVD.size(i,1) = pelvis(i).boundaries.cuboid.SVD.size;             % Size to reference
    pelvisScale.cuboidSVD.scale(i,1) = pelvis(i).boundaries.cuboid.SVD.scale;           % Scaling factor 
    
    % Sphere
    pelvis(i).boundaries.sphere.size = pelvis(i).boundaries.sphere.R / ...
        refSphereR;                                                                     % Size to reference (percent)
    pelvis(i).boundaries.sphere.scale = 1 / pelvis(i).boundaries.sphere.size;           % Scaling factor
    pelvisScale.sphere.R(i,1) = pelvis(i).boundaries.sphere.R;                          % Sphere radius
    pelvisScale.sphere.size(i,1) = pelvis(i).boundaries.sphere.size;                    % Size to reference
    pelvisScale.sphere.scale(i,1) = pelvis(i).boundaries.sphere.scale;                  % Scaling factor 
   
    % Landmarks all
    pelvis(i).boundaries.landmarks.all.size = pelvis(i).boundaries.landmarks.all.distance ./ ...
        refLMdistance;                                                                  % Size to reference (percent), per landmark distance
    pelvis(i).boundaries.landmarks.all.sizeMean = nanmean(pelvis(i).boundaries.landmarks.all.size); % Size to reference (percent)  
    pelvis(i).boundaries.landmarks.all.scaleMean = 1 / pelvis(i).boundaries.landmarks.all.sizeMean; % Scaling factor  
    pelvisScale.LMall.distance(i,:) = pelvis(i).boundaries.landmarks.all.distance;      % Distance
    pelvisScale.LMall.size(i,:) = pelvis(i).boundaries.landmarks.all.size;              % Size to reference, per landmark distance
    pelvisScale.LMall.sizeMean(i,1) = pelvis(i).boundaries.landmarks.all.sizeMean;      % Size to reference
    pelvisScale.LMall.scaleMean(i,1) = pelvis(i).boundaries.landmarks.all.scaleMean;    % Scaling factor
    
    % Landmarks with centre (Acentre)
    pelvis(i).boundaries.landmarks.centre.size = pelvis(i).boundaries.landmarks.centre.distance ./ ...
        refLMdistanceCentre;                                                            % Size to reference (percent), per landmark distance
    pelvis(i).boundaries.landmarks.centre.sizeMean = nanmean(pelvis(i).boundaries.landmarks.centre.size); % Size to reference (percent)
    pelvis(i).boundaries.landmarks.centre.scaleMean = 1 / pelvis(i).boundaries.landmarks.centre.sizeMean; % Scaling factor
    pelvisScale.LMcentre.distance(i,:) = pelvis(i).boundaries.landmarks.centre.distance;% Distance
    pelvisScale.LMcentre.size(i,:) = pelvis(i).boundaries.landmarks.centre.size;        % Size to reference, per landmark distance
    pelvisScale.LMcentre.sizeMean(i,1) = pelvis(i).boundaries.landmarks.centre.sizeMean;% Size to reference
    pelvisScale.LMcentre.scaleMean(i,1) = pelvis(i).boundaries.landmarks.centre.scaleMean;% Scaling factor
    
end

% Additional info to landmarks: mean of LM distance (no scaling)
pelvisScale.LMall.LMdistanceMean = nanmean(pelvisScale.LMall.distance,1); 

% Complete
% Size to reference
pelvisScale.size(:,1) = pelvisScale.LMall.sizeMean;     % ground truth 1
pelvisScale.size(:,2) = pelvisScale.LMcentre.sizeMean;  % ground truth 2
pelvisScale.size(:,3) = pelvisScale.cuboidHull.size;
pelvisScale.size(:,4) = pelvisScale.cuboidSVD.size;
pelvisScale.size(:,5) = pelvisScale.sphere.size;
% Scaling factor
pelvisScale.scale(:,1) = pelvisScale.LMall.scaleMean;   % ground truth 1
pelvisScale.scale(:,2) = pelvisScale.LMcentre.scaleMean;% ground truth 2
pelvisScale.scale(:,3) = pelvisScale.cuboidHull.scale;
pelvisScale.scale(:,4) = pelvisScale.cuboidSVD.scale;
pelvisScale.scale(:,5) = pelvisScale.sphere.scale;

disp('Scaling factor / Size comparison to reference calculated')

%% Deviation/difference of scaling/size methods

% Deviation/difference (absolut) to ground truth 1
for i=1:5
    pelvisScale.sizeDiffabs(:,i) = pelvisScale.size(:,i) - pelvisScale.size(:,1);
    pelvisScale.scaleDiffabs(:,i) = pelvisScale.scale(:,i) - pelvisScale.scale(:,1);
    pelvisScale.meanSizeDiffabs(:,i) = mean(abs(nonzeros(pelvisScale.sizeDiffabs(:,i))));
    pelvisScale.meanScaleDiffabs(:,i) = mean(abs(nonzeros(pelvisScale.scaleDiffabs(:,i))));
    pelvisScale.stdSizeDiffabs(:,i) = std(abs(nonzeros(pelvisScale.sizeDiffabs(:,i))));
    pelvisScale.stdScaleDiffabs(:,i) = std(abs(nonzeros(pelvisScale.scaleDiffabs(:,i))));
end

% Deviation/difference (relativ) to ground truth 1
for i=1:5
    pelvisScale.sizeDiff(:,i) = abs(pelvisScale.sizeDiffabs(:,i)) ./ pelvisScale.size(:,1);
    pelvisScale.sizeDiffs(:,i) = pelvisScale.sizeDiffabs(:,i) ./ pelvisScale.size(:,1);     % with sign
    pelvisScale.scaleDiff(:,i) = abs(pelvisScale.scaleDiffabs(:,i)) ./ pelvisScale.scale(:,1);
    pelvisScale.scaleDiffs(:,i) = pelvisScale.scaleDiffabs(:,i) ./ pelvisScale.scale(:,1);  % with sign
    pelvisScale.meanSizeDiff(:,i) = mean(abs(nonzeros(pelvisScale.sizeDiff(:,i))));
    pelvisScale.meanScaleDiff(:,i) = mean(abs(nonzeros(pelvisScale.scaleDiff(:,i))));
    pelvisScale.stdSizeDiff(:,i) = std(abs(nonzeros(pelvisScale.sizeDiff(:,i))));
    pelvisScale.stdScaleDiff(:,i) = std(abs(nonzeros(pelvisScale.scaleDiff(:,i))));
end


% Deviation/difference (absolut) to ground truth 2
for i=2:5
    pelvisScale.sizeDiff2abs(:,i-1) = pelvisScale.size(:,i) - pelvisScale.size(:,2);
    pelvisScale.scaleDiff2abs(:,i-1) = pelvisScale.scale(:,i) - pelvisScale.scale(:,2);
    pelvisScale.meanSizeDiff2abs(:,i-1) = mean(abs(nonzeros(pelvisScale.sizeDiff2abs(:,i-1))));
    pelvisScale.meanScaleDiff2abs(:,i-1) = mean(abs(nonzeros(pelvisScale.scaleDiff2abs(:,i-1))));
    pelvisScale.stdSizeDiff2abs(:,i-1) = std(abs(nonzeros(pelvisScale.sizeDiff2abs(:,i-1))));
    pelvisScale.stdScaleDiff2abs(:,i-1) = std(abs(nonzeros(pelvisScale.scaleDiff2abs(:,i-1))));
end

% Deviation/difference (relativ) to ground truth 2
for i=2:5
    pelvisScale.sizeDiff2(:,i-1) = abs(pelvisScale.sizeDiff2abs(:,i-1)) ./ pelvisScale.size(:,2);
    pelvisScale.sizeDiff2s(:,i-1) = pelvisScale.sizeDiff2abs(:,i-1) ./ pelvisScale.size(:,2);   % with sign
    pelvisScale.scaleDiff2(:,i-1) = abs(pelvisScale.scaleDiff2abs(:,i-1)) ./ pelvisScale.scale(:,2);
    pelvisScale.scaleDiff2s(:,i-1) = pelvisScale.scaleDiff2abs(:,i-1) ./ pelvisScale.scale(:,2);% with sign
    pelvisScale.meanSizeDiff2(:,i-1) = mean(abs(nonzeros(pelvisScale.sizeDiff2(:,i-1))));
    pelvisScale.meanScaleDiff2(:,i-1) = mean(abs(nonzeros(pelvisScale.scaleDiff2(:,i-1))));
    pelvisScale.stdSizeDiff2(:,i-1) = std(abs(nonzeros(pelvisScale.sizeDiff2(:,i-1))));
    pelvisScale.stdScaleDiff2(:,i-1) = std(abs(nonzeros(pelvisScale.scaleDiff2(:,i-1))));
end

disp('Deviation/difference of scaling/size methods calculated')

%% Display the deviation/difference of size methods

% Size: deviation/difference (relativ) to ground truth 1
figure
bar(1:4,pelvisScale.meanSizeDiff(2:end))
hold on
errorbar(1:4,pelvisScale.meanSizeDiff(2:end),pelvisScale.stdSizeDiff(2:end),...
    pelvisScale.stdSizeDiff(2:end),'LineStyle','none')
% Format and display properties
xticklabels({'LMcentre','CuboidHull','CuboidSVD','Sphere'})
grid on

% Picture BMT paper 2022
figure
bar(2:4,pelvisScale.meanSizeDiff(3:end),'FaceColor',[0 0.396 0.729])
hold on
errorbar(2:4,pelvisScale.meanSizeDiff(3:end),pelvisScale.stdSizeDiff(3:end),...
    pelvisScale.stdSizeDiff(3:end),'LineStyle','none','MarkerFaceColor',[0.89 0.447 0.133],'LineWidth',1.5)
% Format and display properties
xticklabels({'Box Hull','Box SVD','Sphere'})
grid on
xlabel('Scaling Method')
ylabel('Difference in percent')

figure
plot(1:(dataCount-1),pelvisScale.sizeDiffs(2:end,2),'-o') % LMcentre
hold on
plot(1:(dataCount-1),pelvisScale.sizeDiffs(2:end,3),'-o') % Cuboid hull
hold on
plot(1:(dataCount-1),pelvisScale.sizeDiffs(2:end,4),'-o') % Cuboid SVD
hold on
plot(1:(dataCount-1),pelvisScale.sizeDiffs(2:end,5),'-o') % Sphere
hold on
yline(0)
% Format and display properties
grid on
legend('LMcentre','CuboidHull','CuboidSVD','Sphere')
xlabel('Pelvis Number')
ylabel('Size Difference (relative)')

figure
plot(1,pelvisScale.sizeDiffs(2:end,2),'b-o') % LMcentre
hold on
plot(2,pelvisScale.sizeDiffs(2:end,3),'g-o') % Cuboid hull
hold on
plot(3,pelvisScale.sizeDiffs(2:end,4),'r-o') % Cuboid SVD
hold on
plot(4,pelvisScale.sizeDiffs(2:end,5),'c-o') % Sphere
hold on
yline(0)
% Format and display properties
grid on
xlabel('Size comparison method')
ylabel('Size Difference (relative)')


% Size: deviation/difference (relativ) to ground truth 2
figure
bar(1:3,pelvisScale.meanSizeDiff2(2:end))
hold on
errorbar(1:3,pelvisScale.meanSizeDiff2(2:end),pelvisScale.stdSizeDiff2(2:end),...
    pelvisScale.stdSizeDiff2(2:end),'LineStyle','none')
% Format and display properties
xticklabels({'CuboidHull','CuboidSVD','Sphere'})
grid on

figure
plot(1:(dataCount-1),pelvisScale.sizeDiff2s(2:end,2),'-o') % Cuboid hull
hold on
plot(1:(dataCount-1),pelvisScale.sizeDiff2s(2:end,3),'-o') % Cuboid SVD
hold on
plot(1:(dataCount-1),pelvisScale.sizeDiff2s(2:end,4),'-o') % Sphere
hold on
yline(0)
% Format and display properties
grid on
legend('CuboidHull','CuboidSVD','Sphere')
xlabel('Pelvis Number')
ylabel('Size Difference (relative)')

figure
plot(1,pelvisScale.sizeDiff2s(2:end,2),'g-o') % Cuboid hull
hold on
plot(2,pelvisScale.sizeDiff2s(2:end,3),'r-o') % Cuboid SVD
hold on
plot(3,pelvisScale.sizeDiff2s(2:end,4),'c-o') % Sphere
hold on
yline(0)
% Format and display properties
grid on
xlabel('Size comparison method')
ylabel('Size Difference (relative)')


%% Display the deviation/difference of scaling methods

% Scale: deviation/difference (absolut) to ground truth 1
figure
bar(1:4,pelvisScale.meanScaleDiff(2:end))
hold on
errorbar(1:4,pelvisScale.meanScaleDiff(2:end),pelvisScale.stdScaleDiff(2:end),...
    pelvisScale.stdScaleDiff(2:end),'LineStyle','none')
% Format and display properties
xticklabels({'LMcentre','CuboidHull','CuboidSVD','Sphere'})
grid on

figure
plot(1:(dataCount-1),pelvisScale.scaleDiffs(2:end,2),'-o') % LMcentre
hold on
plot(1:(dataCount-1),pelvisScale.scaleDiffs(2:end,3),'-o') % Cuboid hull
hold on
plot(1:(dataCount-1),pelvisScale.scaleDiffs(2:end,4),'-o') % Cuboid SVD
hold on
plot(1:(dataCount-1),pelvisScale.scaleDiffs(2:end,5),'-o') % Sphere
hold on
yline(0)
% Format and display properties
grid on
legend('LMcentre','CuboidHull','CuboidSVD','Sphere')
xlabel('Pelvis Number')
ylabel('Scale Difference (relative)')

figure
plot(1,pelvisScale.scaleDiffs(2:end,2),'b-o') % LMcentre
hold on
plot(2,pelvisScale.scaleDiffs(2:end,3),'g-o') % Cuboid hull
hold on
plot(3,pelvisScale.sizeDiffs(2:end,4),'r-o') % Cuboid SVD
hold on
plot(4,pelvisScale.scaleDiffs(2:end,5),'c-o') % Sphere
hold on
yline(0)
% Format and display properties
grid on
xlabel('Scale comparison method')
ylabel('Scale Difference (relative)')

% Scale: deviation/difference (absolut) to ground truth 2
figure
bar(1:3,pelvisScale.meanScaleDiff2(2:end))
hold on
errorbar(1:3,pelvisScale.meanScaleDiff2(2:end),pelvisScale.stdScaleDiff2(2:end),...
    pelvisScale.stdScaleDiff2(2:end),'LineStyle','none')
% Format and display properties
xticklabels({'CuboidHull','CuboidSVD','Sphere'})
grid on

figure
plot(1:(dataCount-1),pelvisScale.scaleDiff2s(2:end,2),'-o') % Cuboid hull
hold on
plot(1:(dataCount-1),pelvisScale.scaleDiff2s(2:end,3),'-o') % Cuboid SVD
hold on
plot(1:(dataCount-1),pelvisScale.scaleDiff2s(2:end,4),'-o') % Sphere
hold on
yline(0)
% Format and display properties
grid on
legend('CuboidHull','CuboidSVD','Sphere')
xlabel('Pelvis Number')
ylabel('Scale Difference (relative)')

figure
plot(1,pelvisScale.scaleDiff2s(2:end,2),'g-o') % Cuboid hull
hold on
plot(2,pelvisScale.scaleDiff2s(2:end,3),'r-o') % Cuboid SVD
hold on
plot(3,pelvisScale.scaleDiff2s(2:end,4),'c-o') % Sphere
hold on
yline(0)
% Format and display properties
grid on
xlabel('Scale comparison method')
ylabel('Scale Difference (relative)')

%% Dependence of the landmark on the size/scale factor

% Define landmarks and their data
landmarks.names = {'asis', 'tb', 'psis', 'piis', 'si', 'acentre', 'aiis', 'iim', ...
    'iimin', 'spp', 'spd', 'ti', 'fo', 'ci'};
% Landmark distance combis
landmarks.combis = pelvis(1,1).boundaries.landmarks.all.combis;
landmarks.numCombis = size(landmarks.combis,1);


% Pearson Correlation 
% Pearson Correlation Coefficient: For each landmark distance, calculate the Pearson correlation coefficient between
% the landmark distance and the size/scale factor across all pelvis data. A higher absolute value of the coefficient 
% indicates a stronger linear relationship.

% Initialize arrays
landmarks.corr.corrSize = zeros(1, landmarks.numCombis);
landmarks.corr.corrScale = zeros(1, landmarks.numCombis);
landmarks.corr.impactSize = zeros(length(landmarks.names), 1);
landmarks.corr.impactScale = zeros(length(landmarks.names), 1);

% Pearson Correlation
for i = 1:landmarks.numCombis
    % Correlation for size/scale, ignoring NaN values (landmark combination)
    landmarks.corr.corrSize(i) = corr(pelvisScale.LMall.distance(:, i), pelvisScale.LMall.sizeMean, 'rows', 'complete');
    landmarks.corr.corrScale(i) = corr(pelvisScale.LMall.distance(:, i), pelvisScale.LMall.scaleMean, 'rows', 'complete');

    % Accumulating the absolute correlation value for each involved landmark (individual landmark)
    landmarks.corr.impactSize(landmarks.combis(i, :)) = landmarks.corr.impactSize(landmarks.combis(i, :)) + ...
    abs(landmarks.corr.corrSize(i));
    landmarks.corr.impactScale(landmarks.combis(i, :)) = landmarks.corr.impactScale(landmarks.combis(i, :)) + ...
    abs(landmarks.corr.corrScale(i));
end

% Landmark combination/distance
% Sorting landmark combinations based on their correlations
[landmarks.corr.sortedSizeCorr, landmarks.corr.idxSize] = sort(abs(landmarks.corr.corrSize), 'descend');
[landmarks.corr.sortedScaleCorr, landmarks.corr.idxScale] = sort(abs(landmarks.corr.corrScale), 'descend');
% Sorting sorted landmark combinations in arrays
landmarks.corr.sortedSizeArray = cell(landmarks.numCombis, 3);
landmarks.corr.sortedScaleArray = cell(landmarks.numCombis, 3);
for i = 1:landmarks.numCombis
    landmarks.corr.sortedSizeArray{i, 1} = landmarks.names{landmarks.combis(landmarks.corr.idxSize(i), 1)};
    landmarks.corr.sortedSizeArray{i, 2} = landmarks.names{landmarks.combis(landmarks.corr.idxSize(i), 2)};
    landmarks.corr.sortedSizeArray{i, 3} = landmarks.corr.sortedSizeCorr(i);
    landmarks.corr.sortedScaleArray{i, 1} = landmarks.names{landmarks.combis(landmarks.corr.idxScale(i), 1)};
    landmarks.corr.sortedScaleArray{i, 2} = landmarks.names{landmarks.combis(landmarks.corr.idxScale(i), 2)};
    landmarks.corr.sortedScaleArray{i, 3} = landmarks.corr.sortedScaleCorr(i);
end
% Landmarks combination with the highest absolute correlation 
[~, landmarks.corr.maxIndexSize] = max(abs(landmarks.corr.corrSize));
[~, landmarks.corr.maxIndexScale] = max(abs(landmarks.corr.corrScale));
disp(['Landmark combination determining size the most: ', ...
    landmarks.names{landmarks.combis(landmarks.corr.maxIndexSize, 1)}, ' - ', ...
    landmarks.names{landmarks.combis(landmarks.corr.maxIndexSize, 2)}]);
disp(['Landmark combination determining scale the most: ', ...
    landmarks.names{landmarks.combis(landmarks.corr.maxIndexScale, 1)}, ' - ', ...
    landmarks.names{landmarks.combis(landmarks.corr.maxIndexScale, 2)}]);
% auxiliary variables
clear landmarks.corr.idxSize landmarks.corr.idxScale landmarks.corr.maxIndexSize landmarks.corr.maxIndexScale 

% Individual landmark (2 calculations)
% Sorting individual landmarks based on their impact
[~, landmarks.corr.idxSizeImpact] = sort(landmarks.corr.impactSize, 'descend');
[~, landmarks.corr.idxScaleImpact] = sort(landmarks.corr.impactScale, 'descend');
landmarks.corr.sortedLandmarksSizeImpact = landmarks.names(landmarks.corr.idxSizeImpact).';
landmarks.corr.sortedLandmarksScaleImpact = landmarks.names(landmarks.corr.idxScaleImpact).';
% Sorting individual landmarks based on their appearance at the sorted list of landmarks combination
landmarks.corr.scoresSize = zeros(1, numel(landmarks.names));
landmarks.corr.scoresScale = zeros(1, numel(landmarks.names));
% For Size
for idx = 1:size(landmarks.corr.sortedSizeArray, 1)
    landmarks.corr.combis = landmarks.corr.sortedSizeArray(idx, 1:2);
    for l = landmarks.corr.combis
        landmarks.corr.Idx = find(strcmp(landmarks.names, l{1}));
        landmarks.corr.scoresSize(landmarks.corr.Idx) = landmarks.corr.scoresSize(landmarks.corr.Idx) + ...
            (size(landmarks.corr.sortedSizeArray, 1) - idx);
    end
end
% For Scale
for idx = 1:size(landmarks.corr.sortedScaleArray, 1)
    landmarks.corr.combis = landmarks.corr.sortedScaleArray(idx, 1:2);
    for l = landmarks.corr.combis
        landmarks.corr.Idx = find(strcmp(landmarks.names, l{1}));
        landmarks.corr.scoresScale(landmarks.corr.Idx) = landmarks.corr.scoresScale(landmarks.corr.Idx) + ...
            (size(landmarks.corr.sortedScaleArray, 1) - idx);
    end
end
% Sort the scores
[~, landmarks.corr.idxSize] = sort(landmarks.corr.scoresSize, 'descend');
[~, landmarks.corr.idxScale] = sort(landmarks.corr.scoresScale, 'descend');
% Extract the sorted landmarks
landmarks.corr.sortedLandmarksSize = landmarks.names(landmarks.corr.idxSize)';
landmarks.corr.sortedLandmarksScale = landmarks.names(landmarks.corr.idxScale)';
disp(['Landmark determining size the most: ', landmarks.corr.sortedLandmarksSize{1}]);
disp(['Landmark determining scale the most: ', landmarks.corr.sortedLandmarksScale{1}]);
% auxiliary variables
clear landmarks.corr.idxSizeImpact landmarks.corr.idxScaleImpact landmarks.corr.scoresSize landmarks.corr.scoresScale
clear landmarks.corr.combis landmarks.corr.Idx


%% Mean points of the pelves
parfor i=1:dataCount
    pelvisPoints(i,1) = size(pelvis(i).import.processed.vertices,1);
end

pelvisPointsMean = mean(pelvisPoints);

%% Display all size methods

% Size to reference (percent)
figure
bar(pelvisScale.size(2:end,:))
legend('Landmarks','Landmarks Acentre','Cuboid Hull','Cuboid SVD','Sphere')
grid on


%% Dislpay the distribution (histogram) of the pelvis size

% Landmarks
figure
histogram(pelvisScale.LMall.sizeMean(:,1),10)
grid on

figure
edges = (0.7:0.015:1.3);
histogram(pelvisScale.LMall.sizeMean(:,1),edges)
grid on

