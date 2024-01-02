
% @author: Svenja (st100333@stud.uni-stuttgart.de)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a simple function that approximates the minimal Bounding Box of 
% a given 3D dataset.
%
% @Input: 
%           data: nx3 matrix with the x,y,z coordinates of n points
% @Output: 
%           cornerpoints: x,y,z coordinates of the 8 corner points of the
%           oriented bounding box
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cornerpoints] = calc_OriBoundingBox(data)
[n,dim] = size(data);
% Calculate the convexhull of the dataset, the algorithm may work without
% this step but for the most datssets the convexhull gives a better
% approximation to the optimal bounding box
convH = convhull(data(:,1),data(:,2),data(:,3));
% The matlabfunction "convhull" gives a triangulation, it's only necessary
% to know the points on the convex hull.
convH_points = data(convH(:),:);
% ceter the point to calculate the covariance matrix
nK = length(convH_points(:,1));
C = [convH_points(:,1)- sum(convH_points(:,1))/nK,convH_points(:,2)-sum(convH_points(:,2))/nK,convH_points(:,3)-sum(convH_points(:,3))/nK];
cov = C'*C; 
% Use the singular value decomposition to get the principal components of
% data matrix
[U,V,D] = svd(cov);
% Transform the basis to the coordinate system spanned by the eigenvectors
% of V
I = [1 0 0; 0 1 0; 0 0 1];
B_traf = U'\I;
data_traf = data;
for i = 1:n
    data_traf(i,:) = data(i,:)*B_traf;
end
% Calculate the axis assigned bounding box of the rotated data
cornerpoints = zeros(8,3);
x = data_traf(:,1);
y = data_traf(:,2);
z = data_traf(:,3);
% Structure of box (original)
% cornerpoints(1,:) = [min(x), min(y), min(z)];
% cornerpoints(2,:) = [max(x), min(y), min(z)];
% cornerpoints(3,:) = [max(x), max(y), min(z)];
% cornerpoints(4,:) = [min(x), max(y), min(z)];
% cornerpoints(5,:) = [min(x), max(y), max(z)];
% cornerpoints(6,:) = [max(x), max(y), max(z)];
% cornerpoints(7,:) = [max(x), min(y), max(z)];
% cornerpoints(8,:) = [min(x), min(y), max(z)];
% Structure of box (as in minboundbox function)
% 1 xmin,ymin,zmin
% 2 xmax,ymin,zmin
% 3 xmax,ymax,zmin
% 4 xmin,ymax,zmin
% 5 xmin,ymin,zmax
% 6 xmax,ymin,zmax
% 7 xmax,ymax,zmax
% 8 xmin,ymax,zmax
cornerpoints(1,:) = [min(x), min(y), min(z)];
cornerpoints(2,:) = [max(x), min(y), min(z)];
cornerpoints(3,:) = [max(x), max(y), min(z)];
cornerpoints(4,:) = [min(x), max(y), min(z)];
cornerpoints(5,:) = [min(x), min(y), max(z)];
cornerpoints(6,:) = [max(x), min(y), max(z)];
cornerpoints(7,:) = [max(x), max(y), max(z)];
cornerpoints(8,:) = [min(x), max(y), max(z)];
% transform the basis back to the cartesian coordinate system
for j = 1:8
    cornerpoints(j,:) = cornerpoints(j,:)*B_traf';
end
end
