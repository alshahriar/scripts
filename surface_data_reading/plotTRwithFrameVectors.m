function h = plotTRwithFrameVectors(TR, O, ex, ey, ez, varargin)
%PLOTTRWITHFRAMEVECTORS Plot triangulated surface + an existing local frame (O, ex, ey, ez).
%
% h = plotTRwithFrameVectors(TR, O, ex, ey, ez, 'Name',Value,...)
%
% Inputs
%   TR : triangulation
%   O  : 1x3 origin (global)
%   ex,ey,ez : 1x3 unit basis vectors in global coordinates
%
% Name-value options
%   'FaceAlpha' : default 0.35
%   'ShowEdges' : default false
%   'AxisScale' : length of axes arrows in data units (default auto)
%   'LineWidth' : default 2
%   'Title'     : default 'Surface with Local Coordinate Frame'

p = inputParser;
p.addParameter('FaceAlpha', 0.35, @(x) isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
p.addParameter('ShowEdges', false, @(x) islogical(x)&&isscalar(x));
p.addParameter('AxisScale', [], @(x) isempty(x) || (isnumeric(x)&&isscalar(x)&&x>0));
p.addParameter('LineWidth', 2, @(x) isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('Title', 'Surface with Local Coordinate Frame', @(s) ischar(s) || isstring(s));
p.parse(varargin{:});

fa  = p.Results.FaceAlpha;
showEdges = p.Results.ShowEdges;
LW  = p.Results.LineWidth;
ttl = string(p.Results.Title);

F = TR.ConnectivityList;
P = TR.Points;

% Auto axis scale: 15% of bounding box diagonal
if isempty(p.Results.AxisScale)
    diagLen = norm(max(P,[],1) - min(P,[],1));
    if diagLen <= 0, diagLen = 1; end
    L = 0.15 * diagLen;
else
    L = p.Results.AxisScale;
end

figure('Color','w');
ax = axes; hold(ax,'on');

h.surface = patch(ax, 'Faces', F, 'Vertices', P, ...
    'FaceColor', [0.8 0.8 0.8], ...
    'FaceAlpha', fa, ...
    'EdgeColor', ternary(showEdges, [0.2 0.2 0.2], 'none'));

axis(ax,'equal'); view(ax,3); grid(ax,'on');
camlight(ax,'headlight'); lighting(ax,'gouraud');
xlabel(ax,'X'); ylabel(ax,'Y'); zlabel(ax,'Z');

% Axes arrows (X red, Y green, Z blue)
h.xAxis = quiver3(ax, O(1),O(2),O(3), L*ex(1), L*ex(2), L*ex(3), 0, ...
    'LineWidth', LW, 'MaxHeadSize', 0.6, 'Color', [1 0 0]);
h.yAxis = quiver3(ax, O(1),O(2),O(3), L*ey(1), L*ey(2), L*ey(3), 0, ...
    'LineWidth', LW, 'MaxHeadSize', 0.6, 'Color', [0 0.6 0]);
h.zAxis = quiver3(ax, O(1),O(2),O(3), L*ez(1), L*ez(2), L*ez(3), 0, ...
    'LineWidth', LW, 'MaxHeadSize', 0.6, 'Color', [0 0 1]);

% Labels
h.xText = text(ax, O(1)+L*ex(1), O(2)+L*ex(2), O(3)+L*ex(3), '  +X', ...
    'FontWeight','bold','Color',[1 0 0]);
h.yText = text(ax, O(1)+L*ey(1), O(2)+L*ey(2), O(3)+L*ey(3), '  +Y', ...
    'FontWeight','bold','Color',[0 0.6 0]);
h.zText = text(ax, O(1)+L*ez(1), O(2)+L*ez(2), O(3)+L*ez(3), '  +Z', ...
    'FontWeight','bold','Color',[0 0 1]);

% Origin marker
h.origin = plot3(ax, O(1), O(2), O(3), 'ko', 'MarkerFaceColor','k', 'MarkerSize', 5);
h.originText = text(ax, O(1), O(2), O(3), '  O', 'FontWeight','bold', 'Color','k');

title(ax, ttl);
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
