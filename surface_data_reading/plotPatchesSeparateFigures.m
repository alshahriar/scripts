function figHandles = plotPatchesSeparateFigures(TR, faceLabel, varargin)
%PLOTPATCHESSEPARATEFIGURES Plot each surface patch in its own figure.
%
% figHandles = plotPatchesSeparateFigures(TR, faceLabel, 'Name',Value,...)
%
% Inputs
%   TR        : triangulation
%   faceLabel : nFaces x 1 patch IDs
%
% Name-value options
%   'PatchIDs'     : list of patch IDs to plot (default: unique(faceLabel,'stable'))
%   'ShowEdges'    : true/false (default false)
%   'EdgeColor'    : edge color if ShowEdges true (default [0.2 0.2 0.2])
%   'FaceAlpha'    : 0..1 (default 1)
%   'FaceColor'    : face color (default [0.8 0.8 0.8])
%   'LabelTitle'   : true/false (default true)
%   'LinkCamera'   : true/false (default true) link camera across figures
%
% Output
%   figHandles : vector of figure handles

p = inputParser;
p.addRequired('TR', @(x) isa(x,'triangulation'));
p.addRequired('faceLabel', @(x) isnumeric(x) && isvector(x));
p.addParameter('PatchIDs', [], @(x) isempty(x) || (isnumeric(x) && isvector(x)));
p.addParameter('ShowEdges', false, @(x) islogical(x) && isscalar(x));
p.addParameter('EdgeColor', [0.2 0.2 0.2], @(x) (isnumeric(x)&&numel(x)==3) || ischar(x) || isstring(x));
p.addParameter('FaceAlpha', 1, @(x) isnumeric(x) && isscalar(x) && x>=0 && x<=1);
p.addParameter('FaceColor', [0.8 0.8 0.8], @(x) (isnumeric(x)&&numel(x)==3) || ischar(x) || isstring(x));
p.addParameter('LabelTitle', true, @(x) islogical(x) && isscalar(x));
p.addParameter('LinkCamera', true, @(x) islogical(x) && isscalar(x));
p.parse(TR, faceLabel, varargin{:});

F = TR.ConnectivityList;
P = TR.Points;
nF = size(F,1);

faceLabel = faceLabel(:);
if numel(faceLabel) ~= nF
    error('faceLabel must have length equal to number of faces (%d).', nF);
end

if isempty(p.Results.PatchIDs)
    patchIDs = unique(faceLabel, 'stable');
else
    patchIDs = unique(p.Results.PatchIDs(:), 'stable');
end

showEd  = p.Results.ShowEdges;
edgeCol = p.Results.EdgeColor;
aFace   = p.Results.FaceAlpha;
fCol    = p.Results.FaceColor;
doTitle = p.Results.LabelTitle;
doLink  = p.Results.LinkCamera;

figHandles = gobjects(numel(patchIDs), 1);
axHandles  = gobjects(numel(patchIDs), 1);

for k = 1:numel(patchIDs)
    pid = patchIDs(k);
    faceIdx = find(faceLabel == pid);

    if isempty(faceIdx)
        continue;
    end

    Fk = F(faceIdx, :);

    % Re-index vertices for this patch for a compact triangulation
    usedV = unique(Fk(:));
    Pk = P(usedV, :);

    map = zeros(max(usedV), 1);
    map(usedV) = 1:numel(usedV);

    Fk2 = map(Fk);

    figHandles(k) = figure('Name', sprintf('Patch %d', pid), 'Color', 'w');
    axHandles(k) = axes(figHandles(k)); hold(axHandles(k), 'on');

    patch(axHandles(k), 'Faces', Fk2, 'Vertices', Pk, ...
        'FaceColor', fCol, ...
        'FaceAlpha', aFace, ...
        'EdgeColor', ternary(showEd, edgeCol, 'none'));

    axis(axHandles(k), 'equal');
    view(axHandles(k), 3);
    camlight(axHandles(k), 'headlight');
    lighting(axHandles(k), 'gouraud');

    if doTitle
        title(axHandles(k), sprintf('Patch ID: %d   Faces: %d', pid, numel(faceIdx)));
    end
end

% Optionally link camera across all figures/axes (rotate one, others follow)
if doLink
    axHandles = axHandles(isgraphics(axHandles));
    if numel(axHandles) >= 2
        linkprop(axHandles, {'CameraPosition','CameraTarget','CameraUpVector','CameraViewAngle'});
    end
end
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
