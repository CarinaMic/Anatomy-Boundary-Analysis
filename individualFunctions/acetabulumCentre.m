%% Matlab initialisation script

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

%% Data input

% Open file explorer and let user select an STL file
[filename, pathname] = uigetfile('*.stl', 'Select an STL file');

if isequal(filename, 0)
   disp('User selected Cancel');
else
   % Load the STL file using triangulation
   obj = stlread(fullfile(pathname, filename));
   
   % Extract vertices (X) from the object
   X = obj.Points;
   
   % Reduce the number of vertices by removing duplicates
   X = unique(X, 'rows');
end

% Direct import 
% X = pcread('...'); % for point cloud
% X = stlread('...'); % for stl import


%% Sphere Fit (least squared)

% Function sphereFIT
% https://de.mathworks.com/matlabcentral/fileexchange/34129-sphere-fit-least-squared
% this fits a sphere to a collection of data using a closed form for the
% solution (opposed to using an array the size of the data set). 
% Minimizes Sum((x-xc)^2+(y-yc)^2+(z-zc)^2-r^2)^2
% x,y,z are the data, xc,yc,zc are the sphere's center, and r is the radius
% Assumes that points are not in a singular configuration, real numbers, ...
% if you have coplanar data, use a circle fit with svd for determining the
% plane, recommended Circle Fit (Pratt method), by Nikolai Chernov
% http://www.mathworks.com/matlabcentral/fileexchange/22643
% Input:
% X: n x 3 matrix of cartesian data
% Outputs:
% Center: Center of sphere 
% Radius: Radius of sphere
% Author:
% Alan Jennings, University of Dayton

A=[mean(X(:,1).*(X(:,1)-mean(X(:,1)))), ...
    2*mean(X(:,1).*(X(:,2)-mean(X(:,2)))), ...
    2*mean(X(:,1).*(X(:,3)-mean(X(:,3)))); ...
    0, ...
    mean(X(:,2).*(X(:,2)-mean(X(:,2)))), ...
    2*mean(X(:,2).*(X(:,3)-mean(X(:,3)))); ...
    0, ...
    0, ...
    mean(X(:,3).*(X(:,3)-mean(X(:,3))))];
A=A+A.';
B=[mean((X(:,1).^2+X(:,2).^2+X(:,3).^2).*(X(:,1)-mean(X(:,1))));...
    mean((X(:,1).^2+X(:,2).^2+X(:,3).^2).*(X(:,2)-mean(X(:,2))));...
    mean((X(:,1).^2+X(:,2).^2+X(:,3).^2).*(X(:,3)-mean(X(:,3))))];
Center=(A\B).';
Radius=sqrt(mean(sum([X(:,1)-Center(1),X(:,2)-Center(2),X(:,3)-Center(3)].^2,2)));

figure
plot3(X(:,1),X(:,2),X(:,3),'r.')
hold on
daspect([1,1,1]);
[x,y,z] = sphere;
hold on
r = 5;
X2 = x * Radius;
Y2 = y * Radius;
Z2 = z * Radius;
surf(X2+Center(1),Y2+Center(2),Z2+Center(3))

%% Build-in Matlab function pcfitsphere
% https://de.mathworks.com/help/vision/ref/pcfitsphere.html

% pcfitsphere: different answers of function pcfitsphere
% Reason: This behavior is due to the use of MSAC, a version of RANSAC algorithm. 
% RANSAC uses a random process (points) to initialize computation of a mathematical model,
% in this case an essential matrix, and to remove outliers from the data.
% For debugging purpose, set the random seed before the "pcfirsphere" command 
% to get a deterministic result using the following command: rng(0);

% The rng function controls the global stream, which determines how the rand, randi, randn, and 
% randperm functions produce a sequence of random numbers.
% Initializes generator with seed 0
rng(0); 

% Fit a sphere to the point cloud 
maxDistance = 0.001; %%% maximum allowable distance from an inlier point to the sphere (depends on unit of your data)
model = pcfitsphere(pointCloud(X),maxDistance);

% Visualize the point cloud and fitted sphere
figure;
pcshow(X);  % Directly use the matrix X
hold on;
plot(model);
title('Fitted Sphere to Point Cloud');

%% Build-in Matlab function pcfitsphere with iterations
% https://de.mathworks.com/help/vision/ref/pcfitsphere.html

% Set parameters
maxDistance = 0.00001;  % Adjust based on your unit

tic
parfor i = 1:1000 % Iterations
    model = pcfitsphere(pointCloud(X), maxDistance);  % Fit directly using X
    modelCenter(i,:) = model.Center;
    modelRadius(i,1) = model.Radius;
    disp(i)
end
toc

% Calculate mean model
modelMean = sphereModel([mean(modelCenter, 1), mean(modelRadius)]);

% Visualize
figure;
pcshow(X); 
hold on;
plot(modelMean);
title('Mean Fitted Sphere to Point Cloud');