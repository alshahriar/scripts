function S = cfx_readWithReadmatrix(filename, info)
% cfx_readWithReadmatrix
%   Read nodes and faces using readmatrix() using precomputed block info.
%
% Inputs:
%   filename - path to CSV/text file
%   info     - struct returned by cfx_getBlockInfo(filename)
%
% Output:
%   S.info   - same info struct
%   S.nodes  - table of node data
%   S.faces  - numeric matrix of faces (empty if none)

    S = struct();
    S.info  = info;
    S.nodes = table();
    S.faces = [];

    % ---- nodes ----
    if ~isempty(info.node.dataStartLine) && info.node.dataEndLine >= info.node.dataStartLine
        rngN = local_makeRange(info.node.dataStartLine, 1, ...
                               info.node.dataEndLine, info.node.nCols);

        M = readmatrix(filename, ...
            'FileType', 'text', ...
            'Delimiter', ',', ...
            'Range', rngN, ...
            'OutputType', 'double');

        S.nodes = array2table(M, 'VariableNames', info.node.varNames);
    end

    % ---- faces (optional) ----
    if ~isempty(info.face.dataStartLine) && ~isempty(info.face.dataEndLine) && ...
       ~isempty(info.face.nCols) && info.face.nCols > 0 && ...
       info.face.dataEndLine >= info.face.dataStartLine

        rngF = local_makeRange(info.face.dataStartLine, 1, ...
                               info.face.dataEndLine, info.face.nCols);

        S.faces = readmatrix(filename, ...
            'FileType', 'text', ...
            'Delimiter', ',', ...
            'Range', rngF, ...
            'OutputType', 'double');
        S.faces = S.faces + 1;
    end
end

% -------- helpers --------

function rng = local_makeRange(r1, c1, r2, c2)
% Create spreadsheet-style range like 'A12:F200'
    rng = sprintf('%s%d:%s%d', local_colLetters(c1), r1, local_colLetters(c2), r2);
end

function s = local_colLetters(n)
% Convert 1->A, 26->Z, 27->AA, etc.
    s = '';
    while n > 0
        r = mod(n-1, 26);
        s = [char('A' + r) s]; %#ok<AGROW>
        n = floor((n-1)/26);
    end
end
