function faceLabelOut = mergePatchGroups(faceLabelIn, mergeGroups, varargin)
%MERGEPATCHGROUPS Merge multiple groups of patch IDs.
%
% faceLabelOut = mergePatchGroups(faceLabelIn, mergeGroups, 'Renumber', true)
%
% Inputs
%   faceLabelIn  : nFaces x 1 labels
%   mergeGroups  : cell array, each cell is vector of patch IDs to merge
%
% Name-value
%   'Renumber'   : true/false (default true)
%   'Targets'    : optional vector same length as mergeGroups specifying target IDs
%                  (default target for each group = min(group))

p = inputParser;
p.addRequired('faceLabelIn', @(x) isnumeric(x) && isvector(x));
p.addRequired('mergeGroups', @(x) iscell(x) && ~isempty(x));
p.addParameter('Renumber', true, @(x) islogical(x) && isscalar(x));
p.addParameter('Targets', [], @(x) isempty(x) || (isnumeric(x) && isvector(x)));
p.parse(faceLabelIn, mergeGroups, varargin{:});

faceLabel = faceLabelIn(:);
u = unique(faceLabel).';

nG = numel(mergeGroups);
targets = p.Results.Targets;
if isempty(targets)
    targets = zeros(nG,1);
    for k = 1:nG
        g = unique(mergeGroups{k}(:)).';
        if isempty(g), error('mergeGroups{%d} is empty.', k); end
        targets(k) = min(g);
    end
else
    targets = targets(:);
    if numel(targets) ~= nG
        error('Targets must have the same length as mergeGroups.');
    end
end

% Validate and merge sequentially
for k = 1:nG
    g = unique(mergeGroups{k}(:)).';
    missing = setdiff(g, u);
    if ~isempty(missing)
        error('Group %d contains invalid patch IDs: %s', k, mat2str(missing));
    end
    faceLabel(ismember(faceLabel, g)) = targets(k);
end

% Optional renumber to 1..K
if p.Results.Renumber
    [~, ~, faceLabel] = unique(faceLabel, 'stable');
end

faceLabelOut = faceLabel;
end
