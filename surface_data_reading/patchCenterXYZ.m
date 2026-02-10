function [centerMean, centerArea, info] = patchCenterXYZ(TR, faceLabel, patchID)
%PATCHCENTERXYZ Compute the XYZ center of a patch.
%
% Inputs
%   TR        : triangulation object
%   faceLabel : nFaces x 1 patch labels
%   patchID   : scalar patch label to compute center for
%
% Outputs
%   centerMean : 1x3 mean of face centroids (unweighted)
%   centerArea : 1x3 area-weighted centroid (recommended)
%   info       : struct with face indices, areas, face centroids

F = TR.ConnectivityList;
P = TR.Points;

faceLabel = faceLabel(:);
idx = find(faceLabel == patchID);
if isempty(idx)
    error('patchID=%d not found in faceLabel.', patchID);
end

Fi = F(idx,:);

% Face centroids
C = (P(Fi(:,1),:) + P(Fi(:,2),:) + P(Fi(:,3),:)) / 3;

% Mean-of-centroids (simple)
centerMean = mean(C, 1);

% Triangle areas
v1 = P(Fi(:,2),:) - P(Fi(:,1),:);
v2 = P(Fi(:,3),:) - P(Fi(:,1),:);
A  = 0.5 * sqrt(sum(cross(v1, v2, 2).^2, 2));   % nFacesInPatch x 1

% Area-weighted centroid (recommended)
Atot = sum(A);
if Atot == 0
    centerArea = centerMean; % fallback for degenerate patch
else
    centerArea = (A.' * C) / Atot; % 1x3
end

% Return extra info if useful
info.faceIdx     = idx;
info.faceCenters = C;
info.faceAreas   = A;
info.totalArea   = Atot;
end
