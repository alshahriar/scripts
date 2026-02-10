function [faceLabel, faceGroups, boundaryEdges] = splitSurfacesByAngle(TR, angleThreshDeg, varargin)
%SPLITSURFACESBYANGLE Segment a triangulation into surface patches using dihedral angle.
%
% Inputs
%   TR              : triangulation object (triangulation(F,P))
%   angleThreshDeg  : threshold in degrees. Adjacent faces with angle <= threshold are connected.
%
% Name-value options (optional)
%   'NonManifold'   : 'allpairs' (default) | 'ignore'
%                     - allpairs: for edges shared by >2 faces, connect any pair passing threshold
%                     - ignore : do not connect across non-manifold edges
%
% Outputs
%   faceLabel    : nFaces x 1 component label for each face (1..nComponents)
%   faceGroups   : cell array, each cell contains face indices for that component
%   boundaryEdges: Mx2 list of vertex IDs for edges that were "cut" (angle > threshold),
%                  only for edges with >=2 incident faces (manifold and non-manifold cases)
%
% Example
%   [lab, groups] = splitSurfacesByAngle(TR, 20);
%   disp(numel(groups))

% ---- Parse options ----
p = inputParser;
p.addRequired('TR', @(x) isa(x,'triangulation'));
p.addRequired('angleThreshDeg', @(x) isnumeric(x) && isscalar(x) && x>=0 && x<=180);
p.addParameter('NonManifold', 'allpairs', @(s) any(strcmpi(s,{'allpairs','ignore'})));
p.parse(TR, angleThreshDeg, varargin{:});
nonManifoldMode = lower(p.Results.NonManifold);

F = TR.ConnectivityList;   % nF x 3
P = TR.Points;             % nV x 3
nF = size(F,1);

% ---- Face normals (unit) ----
v1 = P(F(:,2),:) - P(F(:,1),:);
v2 = P(F(:,3),:) - P(F(:,1),:);
N  = cross(v1, v2, 2);
Ln = sqrt(sum(N.^2, 2));
Ln(Ln==0) = 1;             % avoid divide-by-zero for degenerate faces
N  = N ./ Ln;

% ---- Build edge -> incident faces map ----
E = [F(:,[1 2]); F(:,[2 3]); F(:,[3 1])];       % 3*nF x 2
E = sort(E, 2);                                  % undirected edges
faceIdx = [(1:nF)'; (1:nF)'; (1:nF)'];          % face for each edge row

[uniqE, ~, ic] = unique(E, 'rows');             % group identical edges
facesPerEdge = accumarray(ic, faceIdx, [], @(x){x});

% ---- Build connectivity graph based on angle threshold ----
s = []; t = [];
boundaryEdges = zeros(0,2);

for k = 1:size(uniqE,1)
    faces = facesPerEdge{k};
    m = numel(faces);

    if m < 2
        continue; % boundary edge (open surface), no neighbor to compare
    end

    if m == 2
        i = faces(1); j = faces(2);
        ang = acosd(max(-1,min(1, dot(N(i,:), N(j,:)))));

        if ang <= angleThreshDeg
            s(end+1,1) = i; %#ok<AGROW>
            t(end+1,1) = j; %#ok<AGROW>
        else
            boundaryEdges(end+1,:) = uniqE(k,:); %#ok<AGROW>
        end

    else
        % Non-manifold edge: shared by >2 faces
        if strcmp(nonManifoldMode, 'ignore')
            boundaryEdges(end+1,:) = uniqE(k,:); %#ok<AGROW>
            continue;
        end

        % allpairs mode: connect any pair of incident faces that pass threshold
        for a = 1:m-1
            for b = a+1:m
                i = faces(a); j = faces(b);
                ang = acosd(max(-1,min(1, dot(N(i,:), N(j,:)))));
                if ang <= angleThreshDeg
                    s(end+1,1) = i; %#ok<AGROW>
                    t(end+1,1) = j; %#ok<AGROW>
                else
                    boundaryEdges(end+1,:) = uniqE(k,:); %#ok<AGROW>
                end
            end
        end
    end
end

% ---- Connected components on face adjacency graph ----
if isempty(s)
    faceLabel = (1:nF).';
else
    G = graph(s, t, [], nF);
    faceLabel = conncomp(G).';  % nF x 1
end

% ---- Pack groups ----
nComp = max(faceLabel);
faceGroups = accumarray(faceLabel, (1:nF).', [nComp 1], @(x){x});
end
