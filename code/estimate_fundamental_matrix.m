% Fundamental Matrix Stencil Code
% CS 4476 / 6476: Computer Vision, Georgia Tech
% Written by Henry Hu

% Returns the camera center matrix for a given projection matrix

% 'Points_a' is nx2 matrix of 2D coordinate of points on Image A
% 'Points_b' is nx2 matrix of 2D coordinate of points on Image B
% 'F_matrix' is 3x3 fundamental matrix

% Try to implement this function as efficiently as possible. It will be
% called repeatly for part III of the project

function [ F_matrix ] = estimate_fundamental_matrix(Points_a,Points_b)

ua = Points_a(:, 1);
va = Points_a(:, 2);
oa = ones(size(ua));

ub = Points_b(:, 1);
vb = Points_b(:, 2);

A = [ua.*ub va.*ub ub ua.*vb va.*vb vb ua va oa];
[~, ~, VA] = svd(A);

FV = reshape(VA(:,9), 3, 3)';

[UF, SF, VF] = svd(FV);
D = diag([SF(1,1) SF(2,2) 0]);
F_matrix = UF * D * VF';

end

