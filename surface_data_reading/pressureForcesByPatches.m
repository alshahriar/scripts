function out = pressureForcesByPatches(S, faceLabel, pressure, varargin)
%PRESSUREFORCESBYPATCHES Compute pressure forces (and torque) for each surface patch.
%
% out = pressureForcesByPatches(S, faceLabel, pressure, 'Name',Value,...)
%
% Added outputs:
%   out.T        : Nt x 3 x Npatch torque vectors
%   out.Tx,Ty,Tz : Nt x Npatch torque components
%
% Extra options:
%   'ComputeTorque'   : true/false (default true)
%   'TorqueReference' : [x0 y0 z0] (default [0 0 0])

% ---- Parse wrapper-specific options (and pass-through the rest) ----
p = inputParser;
p.KeepUnmatched = true;

p.addParameter('PatchIDs', [], @(x) isempty(x) || isnumeric(x) || isstring(x));
p.addParameter('PressureLocation', 'vertex', @(s) any(strcmpi(s, {'vertex','face'})));

p.addParameter('ComputeTorque', true, @(x) islogical(x) && isscalar(x));
p.addParameter('TorqueReference', [0 0 0], @(x) isnumeric(x) && numel(x)==3);

% ForceSign is used by the base function; parse it here so torque matches
p.addParameter('ForceSign', -1, @(x) isnumeric(x) && isscalar(x) && (x==1 || x==-1));

p.parse(varargin{:});

patchIDs     = p.Results.PatchIDs;
pressureLoc  = lower(p.Results.PressureLocation);
doTorque     = p.Results.ComputeTorque;
r0           = reshape(p.Results.TorqueReference, 1, 3);
forceSign    = p.Results.ForceSign;

% Convert Unmatched struct to name-value cell (portable)
un = p.Unmatched;
fn = fieldnames(un);
baseArgs = cell(1, 2*numel(fn));
for i = 1:numel(fn)
    baseArgs{2*i-1} = fn{i};
    baseArgs{2*i}   = un.(fn{i});
end

F  = S.faces;
nF = size(F,1);

faceLabel = faceLabel(:);
if numel(faceLabel) ~= nF
    error('faceLabel must have length equal to number of faces (%d).', nF);
end

if isempty(patchIDs)
    patchIDs = unique(faceLabel, 'stable');
else
    patchIDs = unique(patchIDs(:), 'stable');
end

Nt = size(pressure, 1);

% Consistency checks
switch pressureLoc
    case 'vertex'
        Nnodes = sizeNodesCount(S.nodes);
        if size(pressure,2) ~= Nnodes
            error('For PressureLocation=vertex, pressure must be Nt x Nnodes (%d). Got %d.', ...
                Nnodes, size(pressure,2));
        end
    case 'face'
        if size(pressure,2) ~= nF
            error('For PressureLocation=face, pressure must be Nt x Nfaces (%d). Got %d.', ...
                nF, size(pressure,2));
        end
end

Np = numel(patchIDs);

% ---- Allocate outputs ----
out.PatchIDs = patchIDs;

out.F  = zeros(Nt, 3, Np);
out.Fx = zeros(Nt, Np);
out.Fy = zeros(Nt, Np);
out.Fz = zeros(Nt, Np);

out.T  = zeros(Nt, 3, Np);
out.Tx = zeros(Nt, Np);
out.Ty = zeros(Nt, Np);
out.Tz = zeros(Nt, Np);

out.aux = cell(Np,1);
out.faceIdx = cell(Np,1);
out.nodeIdx = cell(Np,1);

for k = 1:Np
    pid = patchIDs(k);
    faceIdx = find(faceLabel == pid);
    if isempty(faceIdx)
        continue;
    end
    out.faceIdx{k} = faceIdx;

    faces_k = F(faceIdx, :);

    % Compact submesh
    usedNodes = unique(faces_k(:));
    out.nodeIdx{k} = usedNodes;

    map = zeros(sizeNodesCount(S.nodes), 1);
    map(usedNodes) = 1:numel(usedNodes);

    faces_k2 = map(faces_k);

    S_k = struct();
    S_k.faces = faces_k2;
    S_k.nodes = subsetNodes(S.nodes, usedNodes);

    % Pressure subset
    switch pressureLoc
        case 'vertex'
            pressure_k = pressure(:, usedNodes);  % Nt x Nnodes_k
        case 'face'
            pressure_k = pressure(:, faceIdx);    % Nt x Nfaces_k
    end

    % ---- Base force computation (your function) ----
    [Fk, Fxk, Fyk, Fzk, auxk] = pressureForcesTriSurface( ...
        S_k, pressure_k, ...
        'PressureLocation', pressureLoc, ...
        'ForceSign', forceSign, ...
        baseArgs{:});

    out.F(:,:,k) = Fk;
    out.Fx(:,k)  = Fxk;
    out.Fy(:,k)  = Fyk;
    out.Fz(:,k)  = Fzk;
    out.aux{k}   = auxk;

    % ---- Torque computation (optional) ----
    if doTorque
        % Geometry from auxk (already uses your XYZColumns choice)
        a = auxk.a;  % Nfaces_k x 3 area-weighted normals

        % Face vertex indices (fallback to S_k.faces if not present)
        if isfield(auxk, 'i1') && isfield(auxk,'i2') && isfield(auxk,'i3')
            i1 = auxk.i1; i2 = auxk.i2; i3 = auxk.i3;
        else
            i1 = S_k.faces(:,1); i2 = S_k.faces(:,2); i3 = S_k.faces(:,3);
        end

        xyz = auxk.xyz; % Nnodes_k x 3
        C = (xyz(i1,:) + xyz(i2,:) + xyz(i3,:)) / 3;  % Nfaces_k x 3 centroids
        r = C - r0;                                   % lever arms

        % Precompute b_f = r_f x a_f
        b = cross(r, a, 2);                            % Nfaces_k x 3

        % Face pressure over time, Nt x Nfaces_k
        switch pressureLoc
            case 'face'
                pFace = pressure_k; % already face pressure
            case 'vertex'
                % Simple, standard face pressure approximation: average of vertices
                pFace = (pressure_k(:, i1) + pressure_k(:, i2) + pressure_k(:, i3)) / 3;
        end

        Tk = (forceSign .* pFace) * b;  % Nt x 3
        out.T(:,:,k) = Tk;
        out.Tx(:,k)  = Tk(:,1);
        out.Ty(:,k)  = Tk(:,2);
        out.Tz(:,k)  = Tk(:,3);
    end
end
end

% ----------------- helpers -----------------
function n = sizeNodesCount(nodes)
if istable(nodes), n = height(nodes);
else, n = size(nodes,1);
end
end

function nodesSub = subsetNodes(nodes, idx)
if istable(nodes), nodesSub = nodes(idx, :);
else, nodesSub = nodes(idx, :);
end
end
