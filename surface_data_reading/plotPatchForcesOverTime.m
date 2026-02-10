function h = plotPatchForcesOverTime(t, out, varargin)
%PLOTPATCHFORCESOVERTIME Plot patch forces vs time.
%
% h = plotPatchForcesOverTime(t, out, 'Name',Value,...)
%
% Inputs
%   t   : Nt x 1 time vector
%   out : struct returned by pressureForcesByPatches()
%
% Name-value options
%   'Components' : cellstr subset of {'Fx','Fy','Fz','Fmag'} (default all)
%   'PatchIDs'   : patch IDs to plot (default out.PatchIDs)
%   'Legend'     : true/false (default true)
%   'LineWidth'  : default 1.2
%   'SeparateFigures' : true/false (default false)  % one figure per component
%
% Output
%   h : struct of figure/axes handles

p = inputParser;
p.addRequired('t', @(x) isnumeric(x) && isvector(x));
p.addRequired('out', @(x) isstruct(x) && isfield(x,'PatchIDs') && isfield(x,'Fx'));
p.addParameter('Components', {'Fx','Fy','Fz','Fmag'}, @(c) iscell(c) && ~isempty(c));
p.addParameter('PatchIDs', [], @(x) isempty(x) || isnumeric(x) || isstring(x));
p.addParameter('Legend', true, @(x) islogical(x) && isscalar(x));
p.addParameter('LineWidth', 1.2, @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('SeparateFigures', false, @(x) islogical(x) && isscalar(x));
p.parse(t, out, varargin{:});

t = t(:);
Nt = numel(t);

patchIDsAll = out.PatchIDs(:);
if isempty(p.Results.PatchIDs)
    patchIDs = patchIDsAll;
    idxP = 1:numel(patchIDsAll);
else
    patchIDs = unique(p.Results.PatchIDs(:), 'stable');
    [tf, idxP] = ismember(patchIDs, patchIDsAll);
    if any(~tf)
        error('Requested PatchIDs not found in out.PatchIDs: %s', mat2str(patchIDs(~tf).'));
    end
end

Fx = out.Fx(:, idxP);
Fy = out.Fy(:, idxP);
Fz = out.Fz(:, idxP);

if size(Fx,1) ~= Nt
    error('Length of t (%d) must match number of time steps in out.Fx (%d).', Nt, size(Fx,1));
end

Fmag = sqrt(Fx.^2 + Fy.^2 + Fz.^2);

compList = p.Results.Components;
LW = p.Results.LineWidth;
doLegend = p.Results.Legend;
sepFig = p.Results.SeparateFigures;

dataMap = struct('Fx',Fx,'Fy',Fy,'Fz',Fz,'Fmag',Fmag);
yLabelMap = struct('Fx','F_x','Fy','F_y','Fz','F_z','Fmag','|F|');

h = struct();
h.fig = gobjects(0);
h.ax  = gobjects(0);

if sepFig
    for i = 1:numel(compList)
        comp = compList{i};
        if ~isfield(dataMap, comp)
            error('Unknown component "%s". Use Fx, Fy, Fz, Fmag.', comp);
        end

        h.fig(end+1,1) = figure('Color','w');
        h.ax(end+1,1)  = axes(h.fig(end)); hold(h.ax(end),'on');

        plot(h.ax(end), t, dataMap.(comp), 'LineWidth', LW);
        grid(h.ax(end),'on');
        xlabel(h.ax(end),'Time, t');
        ylabel(h.ax(end), yLabelMap.(comp));
        title(h.ax(end), sprintf('%s vs time (all selected patches)', yLabelMap.(comp)));

        if doLegend
            legend(h.ax(end), "Patch " + string(patchIDs(:).'), 'Location','best');
        end
    end
else
    % One figure with stacked layout: one axis per component
    h.fig = figure('Color','w');
    tlo = tiledlayout(h.fig, numel(compList), 1, 'TileSpacing','compact', 'Padding','compact');

    for i = 1:numel(compList)
        comp = compList{i};
        if ~isfield(dataMap, comp)
            error('Unknown component "%s". Use Fx, Fy, Fz, Fmag.', comp);
        end

        h.ax(i,1) = nexttile(tlo); hold(h.ax(i),'on');
        plot(h.ax(i), t, dataMap.(comp), 'LineWidth', LW);
        grid(h.ax(i),'on');
        ylabel(h.ax(i), yLabelMap.(comp));
        title(h.ax(i), sprintf('%s vs time', yLabelMap.(comp)));

        if i == numel(compList)
            xlabel(h.ax(i),'Time, t');
        end

        if doLegend && i == 1
            legend(h.ax(i), "Patch " + string(patchIDs(:).'), 'Location','bestoutside');
        end
    end
end
end
