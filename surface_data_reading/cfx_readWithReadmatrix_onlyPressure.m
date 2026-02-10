function S = cfx_readWithReadmatrix_onlyPressure(filename, info, colNoOfPressure)
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
        rngN = local_makeRange(info.node.dataStartLine, colNoOfPressure, ...
                               info.node.dataEndLine, colNoOfPressure);

        M = readmatrix(filename, ...
            'FileType', 'text', ...
            'Delimiter', ',', ...
            'Range', rngN, ...
            'OutputType', 'double');
        S.nodes = array2table(M, 'VariableNames', info.node.varNames(colNoOfPressure));
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
