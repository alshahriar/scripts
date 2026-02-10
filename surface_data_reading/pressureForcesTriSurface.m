function [F, Fx, Fy, Fz, aux] = pressureForcesTriSurface(S, pressure, varargin)
%pressureForcesTriSurface Net pressure force on a triangulated surface vs time.
%
%   [F, Fx, Fy, Fz, aux] = pressureForcesTriSurface(S, pressure)
%   computes the net force vector due to pressure loading on the surface:
%       F(t) = -∫_S p(x,t) n dA  ≈  -Σ_f p_f(t) (n_f A_f)
%
% INPUTS
%   S        : struct with fields
%              - nodes : (Nnodes x >=3) table or numeric array. Must contain XYZ.
%              - faces : (Nfaces x 3) double, 1-based vertex indices.
%   pressure : (Nt x Nnodes) vertex pressure per time step, OR (Nt x Nfaces)
%              face pressure per time step (see 'PressureLocation').
%
% NAME-VALUE OPTIONS
%   'PressureLocation' : 'vertex' (default) or 'face'
%   'XYZColumns'       : vector of 3 column indices into nodes for [x y z]
%                        default: if table, tries x/y/z by name else [1 2 3]
%   'ChunkSize'        : positive integer, default 50
%   'ForceSign'        : -1 (default) or +1 (flip sign if normals are opposite)
%
% OUTPUTS
%   F    : (Nt x 3) [Fx Fy Fz]
%   Fx,Fy,Fz : (Nt x 1) components
%   aux  : struct with geometry precomputations:
%          aux.a          (Nfaces x 3) area-weighted normals (n*A)
%          aux.i1,i2,i3   face vertex indices
%          aux.xyz        (Nnodes x 3) xyz used
%
% NOTES
%   - Assumes S.faces are triangles.
%   - If pressure is in Pa and xyz in m, force is in N.
%
% Example:
%   [F,Fx,Fy,Fz] = pressureForcesTriSurface(S, pressure, ...
%       'PressureLocation','vertex','ChunkSize',100);

% -------------------- Parse options --------------------
p = inputParser;
p.addParameter('PressureLocation','vertex', @(s)ischar(s) || isstring(s));
p.addParameter('XYZColumns',[], @(v)isnumeric(v) && (isempty(v) || numel(v)==3));
p.addParameter('ChunkSize',20, @(v)isnumeric(v) && isscalar(v) && v>=1);
p.addParameter('ForceSign',-1, @(v)isnumeric(v) && isscalar(v) && (v==1 || v==-1));
p.parse(varargin{:});

pressureLocation = lower(string(p.Results.PressureLocation));
XYZColumns       = p.Results.XYZColumns;
chunkSize        = p.Results.ChunkSize;
forceSign        = p.Results.ForceSign;

% -------------------- Extract geometry --------------------
faces = S.faces;
if size(faces,2) ~= 3
    error('S.faces must be Nfaces x 3 (triangles).');
end

% Extract xyz
if istable(S.nodes)
    vn = S.nodes.Properties.VariableNames;
    if ~isempty(XYZColumns)
        xyz = S.nodes{:, XYZColumns};
    else
        % Try variable names x,y,z (case-insensitive), else first 3 columns
        ix = find(strcmpi(vn,'x'),1);
        iy = find(strcmpi(vn,'y'),1);
        iz = find(strcmpi(vn,'z'),1);
        if ~isempty(ix) && ~isempty(iy) && ~isempty(iz)
            xyz = S.nodes{:, [ix iy iz]};
        else
            xyz = S.nodes{:, 1:3};
        end
    end
else
    if ~isempty(XYZColumns)
        xyz = S.nodes(:, XYZColumns);
    else
        xyz = S.nodes(:, 1:3);
    end
end

% Face vertex indices
i1 = faces(:,1);
i2 = faces(:,2);
i3 = faces(:,3);

% Area-weighted normals: a_f = n*A (Nfaces x 3)
v1 = xyz(i1,:);
v2 = xyz(i2,:);
v3 = xyz(i3,:);
a  = 0.5 * cross(v2 - v1, v3 - v1, 2);

% -------------------- Validate pressure shape --------------------
Nt = size(pressure,1);
Nfaces = size(faces,1);
Nnodes = size(xyz,1);

switch pressureLocation
    case "vertex"
        if size(pressure,2) ~= Nnodes
            error('PressureLocation="vertex" requires pressure to be Nt x Nnodes (%d). Got %d.', ...
                Nnodes, size(pressure,2));
        end
    case "face"
        if size(pressure,2) ~= Nfaces
            error('PressureLocation="face" requires pressure to be Nt x Nfaces (%d). Got %d.', ...
                Nfaces, size(pressure,2));
        end
    otherwise
        error('PressureLocation must be "vertex" or "face".');
end

% -------------------- Compute forces (chunked) --------------------
F = zeros(Nt, 3);

for t0 = 1:chunkSize:Nt
    t1 = min(t0 + chunkSize - 1, Nt);

    if pressureLocation == "vertex"
        % p_f = average vertex pressures per face
        pface = (pressure(t0:t1, i1) + pressure(t0:t1, i2) + pressure(t0:t1, i3)) / 3;
    else
        % already face pressure
        pface = pressure(t0:t1, :);
    end

    % F = sign * (pface * a) where default sign = -1
    F(t0:t1,:) = forceSign * (pface * a);
end

Fx = F(:,1); Fy = F(:,2); Fz = F(:,3);

% -------------------- Aux outputs --------------------
aux = struct();
aux.a   = a;
aux.i1  = i1; aux.i2 = i2; aux.i3 = i3;
aux.xyz = xyz;
aux.Nt = Nt; aux.Nnodes = Nnodes; aux.Nfaces = Nfaces;
aux.pressureLocation = char(pressureLocation);
aux.forceSign = forceSign;

end
