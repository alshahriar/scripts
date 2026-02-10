function [faceLabelOut, newIDs, splitInfo] = splitPatchByRadiusXZ(TR, faceLabelIn, patchID, rSplit, varargin)
%SPLITPATCHBYRADIUSXZ Split a patch into inside/outside based on radius in XZ plane.
%
% [faceLabelOut, newIDs, splitInfo] = splitPatchByRadiusXZ(TR, faceLabelIn, patchID, rSplit, ...)
%
% Inputs
%   TR          : triangulation object
%   faceLabelIn : nFaces x 1 patch IDs
%   patchID     : scalar patch ID to split
%   rSplit      : scalar radius threshold in XZ plane
%
% Name-value options
%   'CenterXZ'  : [x0 z0] (default [0 0])
%   'Method'    : 'centroid' (default) | 'vertexFraction'
%   'Fraction'  : used only for 'vertexFraction' (default 0.5)
%                 - face goes "inside" if (insideVertices/3) >= Fraction
%   'InsideID'  : label to assign to inside faces (default: patchID)
%   'OutsideID' : label to assign to outside faces (default: new max+1)
%   'Renumber'  : true/false (default false)
%
% Outputs
%   faceLabelOut : updated labels
%   newIDs       : struct with fields .InsideID, .OutsideID
%   splitInfo    : struct with counts and indices

p = inputParser;
p.addRequired('TR', @(x) isa(x,'triangulation'));
p.addRequired('faceLabelIn', @(x) isnumeric(x) && isvector(x));
p.addRequired('patchID', @(x) isnumeric(x) && isscalar(x));
p.addRequired('rSplit', @(x) isnumeric(x) && isscalar(x) && x>=0);

p.addParameter('CenterXZ', [0 0], @(x) isnumeric(x) && numel(x)==2);
p.addParameter('Method', 'centroid', @(s) any(strcmpi(s, {'centroid','vertexfraction'})));
p.addParameter('Fraction', 0.5, @(x) isnumeric(x) && isscalar(x) && x>=0 && x<=1);
p.addParameter('InsideID', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
p.addParameter('OutsideID', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
p.addParameter('Renumber', false, @(x) islogical(x) && isscalar(x));
p.parse(TR, faceLabelIn, patchID, rSplit, varargin{:});

F = TR.ConnectivityList;
P = TR.Points;
nF = size(F,1);

faceLabel = faceLabelIn(:);
if numel(faceLabel) ~= nF
    error('faceLabelIn must have length equal to number of faces (%d).', nF);
end

if ~any(faceLabel == patchID)
    error('patchID=%d does not exist in faceLabelIn.', patchID);
end

x0 = p.Results.CenterXZ(1);
z0 = p.Results.CenterXZ(2);

method = lower(p.Results.Method);

% Determine inside/outside classification for faces in patchID
inPatch = (faceLabel == patchID);
idxFaces = find(inPatch);

switch method
    case 'centroid'
        C = (P(F(:,1),:) + P(F(:,2),:) + P(F(:,3),:)) / 3; % centroids
        r = sqrt((C(:,1)-x0).^2 + (C(:,3)-z0).^2);
        inside = (r <= rSplit);

    case 'vertexfraction'
        % classify a face based on how many of its vertices are within rSplit
        vx = P(:,1); vz = P(:,3);
        rv = sqrt((vx-x0).^2 + (vz-z0).^2);
        insideV = rv <= rSplit; % per-vertex
        fracInside = (insideV(F(:,1)) + insideV(F(:,2)) + insideV(F(:,3))) / 3;
        inside = (fracInside >= p.Results.Fraction);
end

insidePatch = inPatch & inside;
outsidePatch = inPatch & ~inside;

% Choose IDs
if isempty(p.Results.InsideID)
    insideID = patchID;           % keep inside as original patch by default
else
    insideID = p.Results.InsideID;
end

if isempty(p.Results.OutsideID)
    outsideID = max(faceLabel) + 1;
else
    outsideID = p.Results.OutsideID;
end

% Apply relabeling
faceLabel(insidePatch)  = insideID;
faceLabel(outsidePatch) = outsideID;

% Optional renumber to 1..K (stable)
if p.Results.Renumber
    [~,~,faceLabel] = unique(faceLabel, 'stable');
end

faceLabelOut = faceLabel;

newIDs.InsideID  = insideID;
newIDs.OutsideID = outsideID;

splitInfo.patchID = patchID;
splitInfo.rSplit  = rSplit;
splitInfo.centerXZ = p.Results.CenterXZ;
splitInfo.method  = method;
splitInfo.nInside = nnz(insidePatch);
splitInfo.nOutside = nnz(outsidePatch);
splitInfo.faceIdxInside = find(insidePatch);
splitInfo.faceIdxOutside = find(outsidePatch);
splitInfo.faceIdxOriginal = idxFaces;
end
