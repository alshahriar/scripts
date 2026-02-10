function gifFile = makePressureGifTwoViews(TR, pressure, time, mean_P, directory, varargin)
%makePressureGifTwoViews Create a 2-panel (side-by-side) GIF animation of pressure.
%
%   gifFile = makePressureGifTwoViews(TR, pressure, time, mean_P, directory)
%
% REQUIRED INPUTS
%   TR        : triangulation object (or compatible input to trisurf)
%   pressure  : Nt x Nnodes (per-vertex) OR Nt x Nfaces (per-face) array
%   time      : Nt x 1 (or 1 x Nt) vector
%   mean_P    : scalar used for clim = mean_P ± 5%
%   directory : output folder for the gif
%
% NAME-VALUE OPTIONS
%   'GifName'        : default "animation2.gif"
%   'StartIndex'     : default 400
%   'EndIndex'       : default size(pressure,1)
%   'DelayTime'      : default 0.05
%   'LoopCount'      : default inf
%   'FrameSkip'      : default 1
%   'ViewLeft'       : default [0 -90]
%   'ViewRight'      : default [-180 90]
%   'AxisOff'        : default true
%   'TileSpacing'    : default 'compact'
%   'Padding'        : default 'compact'
%   'CamPosLeft'     : [] (no override) or 1x3 vector
%   'CamPosRight'    : [] (no override) or 1x3 vector
%   'FigSize'        : [] (no override) or [width height] in pixels
%
% OUTPUT
%   gifFile : full path of the written GIF

% -------------------- parse inputs --------------------
p = inputParser;
p.addParameter('GifName', 'animation2.gif', @(s)ischar(s) || isstring(s));
p.addParameter('StartIndex', 400, @(v)isnumeric(v) && isscalar(v) && v>=1);
p.addParameter('EndIndex', [], @(v)isnumeric(v) && isscalar(v) && v>=1);
p.addParameter('DelayTime', 0.05, @(v)isnumeric(v) && isscalar(v) && v>=0);
p.addParameter('LoopCount', inf, @(v)isnumeric(v) && isscalar(v));
p.addParameter('FrameSkip', 1, @(v)isnumeric(v) && isscalar(v) && v>=1);

p.addParameter('ViewLeft',  [0 -90],    @(v)isnumeric(v) && numel(v)==2);
p.addParameter('ViewRight', [-180 90],  @(v)isnumeric(v) && numel(v)==2);
p.addParameter('AxisOff', true, @(v)islogical(v) && isscalar(v));

p.addParameter('TileSpacing', 'compact', @(s)ischar(s) || isstring(s));
p.addParameter('Padding',     'compact', @(s)ischar(s) || isstring(s));

p.addParameter('CamPosLeft',  [], @(v)isnumeric(v) && (isempty(v) || numel(v)==3));
p.addParameter('CamPosRight', [], @(v)isnumeric(v) && (isempty(v) || numel(v)==3));

p.addParameter('FigSize', [], @(v)isnumeric(v) && (isempty(v) || (numel(v)==2 && all(v>0))));
p.parse(varargin{:});

gifName     = string(p.Results.GifName);
iStart      = p.Results.StartIndex;
iEnd        = p.Results.EndIndex;
delayTime   = p.Results.DelayTime;
loopCount   = p.Results.LoopCount;
frameSkip   = p.Results.FrameSkip;

viewLeft    = p.Results.ViewLeft;
viewRight   = p.Results.ViewRight;
axisOff     = p.Results.AxisOff;

tileSpacing = string(p.Results.TileSpacing);
padding     = string(p.Results.Padding);

camPosLeft  = p.Results.CamPosLeft;
camPosRight = p.Results.CamPosRight;
figSize     = p.Results.FigSize;

Nt = size(pressure,1);
if isempty(iEnd); iEnd = Nt; end
iStart = max(1, min(iStart, Nt));
iEnd   = max(1, min(iEnd,   Nt));

time = time(:);
if numel(time) ~= Nt
    error('time must have Nt elements (Nt = size(pressure,1) = %d).', Nt);
end

% -------------------- figure & axes --------------------
fig = figure('Color','w');
tlo = tiledlayout(fig, 1, 2, 'TileSpacing', char(tileSpacing), 'Padding', char(padding));

if ~isempty(figSize)
    % Keep position origin, set size deterministically
    pos = fig.Position;
    fig.Position = [pos(1) pos(2) figSize(1) figSize(2)];
end

% Left axes
ax1 = nexttile(tlo, 1);
h1  = trisurf(TR, pressure(1,:), "edgeColor","none", "Parent", ax1);
axis(ax1,'equal'); axis(ax1,'vis3d');
view(ax1, viewLeft);
lighting(ax1,'phong'); colorbar(ax1);
clim(ax1,[mean_P-0.05*mean_P, mean_P+0.05*mean_P]);
colormap(ax1,"jet");
light(ax1,'Position', [0 0 0],  'Style', 'local', 'Color', 'w');
light(ax1,'Position', [5 5 10], 'Style', 'local', 'Color', 'w');
if ~isempty(camPosLeft); campos(ax1, camPosLeft); end
if axisOff; axis(ax1,'off'); end

% Right axes
ax2 = nexttile(tlo, 2);
h2  = trisurf(TR, pressure(1,:), "edgeColor","none", "Parent", ax2);
axis(ax2,'equal'); axis(ax2,'vis3d');
view(ax2, viewRight);
lighting(ax2,'phong'); colorbar(ax2);
clim(ax2,[mean_P-0.05*mean_P, mean_P+0.05*mean_P]);
colormap(ax2,"jet");
light(ax2,'Position', [0 0 0],  'Style', 'local', 'Color', 'w');
light(ax2,'Position', [5 5 10], 'Style', 'local', 'Color', 'w');
if ~isempty(camPosRight); campos(ax2, camPosRight); end
if axisOff; axis(ax2,'off'); end

% -------------------- GIF setup --------------------
gifFile = fullfile(directory, gifName);
if exist(gifFile,'file'); delete(gifFile); end
writtenFirst = false;

% -------------------- main loop --------------------
for iFile = iStart:iEnd
    pnow = pressure(iFile,:);

    set(h1,'CData',pnow);
    set(h2,'CData',pnow);

    % Keep views fixed
    view(ax1, viewLeft);
    view(ax2, viewRight);

    title(ax1, "time = " + string(time(iFile)) + " [s]");
    title(ax2, "time = " + string(time(iFile)) + " [s]");

    drawnow

    if mod(iFile, frameSkip) ~= 0
        continue
    end

    fr = getframe(fig);
    [im, cm] = rgb2ind(fr.cdata, 256);

    if ~writtenFirst
        imwrite(im, cm, gifFile, 'gif', 'LoopCount', loopCount, 'DelayTime', delayTime);
        writtenFirst = true;
    else
        imwrite(im, cm, gifFile, 'gif', 'WriteMode', 'append', 'DelayTime', delayTime);
    end
end
end
