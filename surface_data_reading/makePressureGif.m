function gifFile = makePressureGif(TR, pressure, time, mean_P, directory, varargin)
%makePressureGif Create a GIF animation of pressure on a triangulated surface.
%
%   gifFile = makePressureGif(TR, pressure, time, mean_P, directory)
%
% REQUIRED INPUTS
%   TR        : triangulation object (or compatible input to trisurf)
%   pressure  : Nt x Nnodes (per-vertex) OR Nt x Nfaces (per-face) array
%   time      : Nt x 1 (or 1 x Nt) vector
%   mean_P    : scalar used for clim = mean_P ± 5%
%   directory : output folder for the gif
%
% NAME-VALUE OPTIONS
%   'GifName'       : default "animation3.gif"
%   'StartIndex'    : default 400
%   'EndIndex'      : default size(pressure,1)
%   'DelayTime'     : default 0.05
%   'LoopCount'     : default inf
%   'FrameSkip'     : default 1
%   'View'          : default [0 180]
%   'AxisOff'       : default true
%   'CamPos'        : default [0.5129 -0.1886 -0.1439]
%   'UseCameraSpin' : default false (uses Rodrigues rotation like your code)
%   'SpinAxis'      : default [1 2 3]
%   'SpinDegPerFrame': default 0.3
%
% OUTPUT
%   gifFile : full path of the written GIF

% -------------------- parse inputs --------------------
p = inputParser;
p.addParameter('GifName', 'animation3.gif', @(s)ischar(s) || isstring(s));
p.addParameter('StartIndex', 400, @(v)isnumeric(v) && isscalar(v) && v>=1);
p.addParameter('EndIndex', [], @(v)isnumeric(v) && isscalar(v) && v>=1);
p.addParameter('DelayTime', 0.05, @(v)isnumeric(v) && isscalar(v) && v>=0);
p.addParameter('LoopCount', inf, @(v)isnumeric(v) && isscalar(v));
p.addParameter('FrameSkip', 1, @(v)isnumeric(v) && isscalar(v) && v>=1);
p.addParameter('View', [0 180], @(v)isnumeric(v) && numel(v)==2);
p.addParameter('AxisOff', true, @(v)islogical(v) && isscalar(v));
p.addParameter('CamPos', [0.5129 -0.1886 -0.1439], @(v)isnumeric(v) && numel(v)==3);
p.addParameter('UseCameraSpin', false, @(v)islogical(v) && isscalar(v));
p.addParameter('SpinAxis', [1 2 3], @(v)isnumeric(v) && numel(v)==3);
p.addParameter('SpinDegPerFrame', 0.3, @(v)isnumeric(v) && isscalar(v));
p.parse(varargin{:});

gifName        = string(p.Results.GifName);
iStart         = p.Results.StartIndex;
iEnd           = p.Results.EndIndex;
delayTime      = p.Results.DelayTime;
loopCount      = p.Results.LoopCount;
frameSkip      = p.Results.FrameSkip;
viewAng        = p.Results.View;
axisOff        = p.Results.AxisOff;
camPos         = p.Results.CamPos;
useSpin        = p.Results.UseCameraSpin;
spinAxis       = p.Results.SpinAxis(:).';
spinDegPerFrame= p.Results.SpinDegPerFrame;

Nt = size(pressure,1);
if isempty(iEnd); iEnd = Nt; end
iStart = max(1, min(iStart, Nt));
iEnd   = max(1, min(iEnd,   Nt));

time = time(:); % column
if numel(time) ~= Nt
    error('time must have Nt elements (Nt = size(pressure,1) = %d).', Nt);
end

% -------------------- initial figure setup --------------------
fig = figure('Color','w');
ax  = axes(fig);

h = trisurf(TR, pressure(1,:), "edgeColor","none", "Parent", ax);
axis(ax,'equal'); axis(ax,'vis3d');

light(ax, 'Position', [0 0 0],  'Style', 'local', 'Color', 'w');
light(ax, 'Position', [5 5 10], 'Style', 'local', 'Color', 'w');

lighting(ax,'phong');
colorbar(ax);
clim(ax,[mean_P-0.05*mean_P, mean_P+0.05*mean_P]);
colormap(ax,"jet");

campos(ax, camPos);

% Spin rotation (optional)
spinAxis = spinAxis / norm(spinAxis);
R = rodriguesRot(spinAxis, deg2rad(spinDegPerFrame));  %#ok<NASGU>

% -------------------- GIF setup --------------------
gifFile = fullfile(directory, gifName);
if exist(gifFile,'file'); delete(gifFile); end
writtenFirst = false;

% -------------------- main loop --------------------
for iFile = iStart:iEnd
    set(h, 'CData', pressure(iFile,:));
    title(ax, "time = " + string(time(iFile)) + " [s]");

    if useSpin
        tgt = camtarget(ax);
        pos = campos(ax);
        up  = camup(ax);
        campos(ax, (R*(pos - tgt).').'+tgt);
        camup(ax,  (R*up.').');
    else
        view(ax, viewAng);
        if axisOff
            axis(ax,'off');
        end
    end

    drawnow;

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
