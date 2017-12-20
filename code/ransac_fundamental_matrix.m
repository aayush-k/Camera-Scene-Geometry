% RANSAC Stencil Code
% CS 4476 / 6476: Computer Vision, Georgia Tech
% Written by Henry Hu

% Find the best fundamental matrix using RANSAC on potentially matching
% points

% 'matches_a' and 'matches_b' are the Nx2 coordinates of the possibly
% matching points from pic_a and pic_b. Each row is a correspondence (e.g.
% row 42 of matches_a is a point that corresponds to row 42 of matches_b.

% 'Best_Fmatrix' is the 3x3 fundamental matrix
% 'inliers_a' and 'inliers_b' are the Mx2 corresponding points (some subset
% of 'matches_a' and 'matches_b') that are inliers with respect to
% Best_Fmatrix.

% For this section, use RANSAC to find the best fundamental matrix by
% randomly sample interest points. You would reuse
% estimate_fundamental_matrix() from part 2 of this assignment.

% If you are trying to produce an uncluttered visualization of epipolar
% lines, you may want to return no more than 30 points for either left or
% right images.

function [ Best_Fmatrix, inliers_a, inliers_b] = ransac_fundamental_matrix(matches_a, matches_b)

matchDim = size(matches_a);
iterations = 1000;
sampleSize = 8;
threshold = 0.02;
maxInliers = 0;
bestInlierInd = 0;
Best_Fmatrix = 0;

for i = 1:iterations

    sample = randsample(matchDim(1), sampleSize);

    F = estimate_fundamental_matrix(matches_a(sample, :), matches_b(sample, :));

    oneRow = ones(1, matchDim(1));
    F_a = F * [matches_a'; oneRow];
    F_b = F * [matches_b'; oneRow];

    F2 = zeros(1, matchDim(1));
    for j = 1:matchDim(1)
        % homogenize
        mB_hom = [matches_b(j,:) 1];
        mA_hom = [matches_a(j,:) 1];
        F2(j) = mB_hom * F * mA_hom';
    end

    % Get Sampson Error and threshold to find inliers
    error = (F2 .^ 2) ./ (F_a(1,:) .^ 2 + F_a(2,:) .^ 2 + F_b(1,:) .^ 2 + F_b(2,:) .^ 2);
    inlierInd = find(abs(error) < threshold);
    numInliers = length(inlierInd);

    % see if better inliers are found
    if (numInliers > maxInliers)
        % better model found
        maxInliers = numInliers;
        bestInlierInd = inlierInd;
        Best_Fmatrix = F;

    end


end

inliers_a = matches_a(bestInlierInd, :);
inliers_b = matches_b(bestInlierInd, :);

end

