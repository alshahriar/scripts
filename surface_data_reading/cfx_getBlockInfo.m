function info = cfx_getBlockInfo(filename)
% cfx_getBlockInfo
%   Detects line ranges and column counts for:
%     - Node data block under [Data]
%     - Face connectivity block under [Faces] (optional)
%
% Returns struct 'info' with fields:
%   info.node.headerLine
%   info.node.dataStartLine
%   info.node.dataEndLine
%   info.node.nCols
%   info.node.headerText
%   info.node.varNames
%
%   info.face.dataStartLine
%   info.face.dataEndLine
%   info.face.nCols
%
% Line numbers are 1-based (MATLAB style).

    if ~isfile(filename)
        error('cfx_getBlockInfo:FileNotFound', 'File not found: %s', filename);
    end

    txt   = fileread(filename);
    lines = regexp(txt, '\r\n|\n|\r', 'split');
    lines = lines(:);
    nL    = numel(lines);

    % Helper lambdas
    trimL = @(s) strtrim(s);
    isTag = @(s) ~isempty(s) && startsWith(trimL(s), '[') && endsWith(trimL(s), ']');

    idxData  = find(strcmpi(trimL(lines), '[Data]'), 1, 'first');
    idxFaces = find(strcmpi(trimL(lines), '[Faces]'), 1, 'first');

    if isempty(idxData)
        error('cfx_getBlockInfo:NoDataTag', 'No [Data] tag found in file.');
    end

    % ---- NODE header line: first non-empty line after [Data]
    headerLine = local_nextNonEmpty(lines, idxData + 1);
    if isempty(headerLine)
        error('cfx_getBlockInfo:MissingHeader', 'Found [Data] but no header line after it.');
    end

    headerText = trimL(lines{headerLine});
    varNames   = local_parseHeader(headerText);
    nColsNode  = numel(varNames);

    % ---- NODE data start: first non-empty (and not a tag) after header line
    dataStart = local_nextNonEmpty(lines, headerLine + 1);
    while ~isempty(dataStart) && isTag(lines{dataStart})
        dataStart = local_nextNonEmpty(lines, dataStart + 1);
    end

    % ---- NODE data end: just before [Faces] if present, else before next tag or EOF
    if ~isempty(idxFaces) && idxFaces > headerLine
        dataEnd = idxFaces - 1;
    else
        dataEnd = nL;
    end
    dataEnd = local_lastNonEmptyBeforeTag(lines, dataStart, dataEnd);

    % Pack node info
    info = struct();
    info.node = struct();
    info.node.headerLine    = headerLine;
    info.node.dataStartLine = dataStart;
    info.node.dataEndLine   = dataEnd;
    info.node.nCols         = nColsNode;
    info.node.headerText    = headerText;
    info.node.varNames      = lower(varNames);

    % ---- FACE block (optional)
    info.face = struct();
    info.face.dataStartLine = [];
    info.face.dataEndLine   = [];
    info.face.nCols         = [];

    if ~isempty(idxFaces)
        fStart = local_nextNonEmpty(lines, idxFaces + 1);
        while ~isempty(fStart) && isTag(lines{fStart})
            fStart = local_nextNonEmpty(lines, fStart + 1);
        end

        fEnd = nL;
        fEnd = local_lastNonEmptyBeforeTag(lines, fStart, fEnd);

        % Determine face columns from first numeric line in faces
        firstFace = local_nextNumericLine(lines, fStart, fEnd);
        if isempty(firstFace)
            % Faces tag exists but no numeric face data found
            info.face.dataStartLine = fStart;
            info.face.dataEndLine   = fEnd;
            info.face.nCols         = 0;
        else
            t = trimL(lines{firstFace});
            info.face.dataStartLine = fStart;
            info.face.dataEndLine   = fEnd;
            info.face.nCols         = count(t, ',') + 1;
        end
    end
end

% ---------- helpers ----------

function k = local_nextNonEmpty(lines, startIdx)
    k = [];
    for i = startIdx:numel(lines)
        if ~isempty(strtrim(lines{i}))
            k = i;
            return;
        end
    end
end

function k = local_nextNumericLine(lines, startIdx, endIdx)
% Numeric line: begins with digit, +, -, or .
    k = [];
    if isempty(startIdx), return; end
    endIdx = min(endIdx, numel(lines));
    for i = startIdx:endIdx
        s = strtrim(lines{i});
        if isempty(s), continue; end
        if ~isempty(regexp(s, '^[\+\-\.0-9]', 'once'))
            k = i; return;
        end
    end
end

function endLine = local_lastNonEmptyBeforeTag(lines, startLine, endLine)
% Shrinks endLine upward over blanks; also stops earlier if a tag line appears.
    if isempty(startLine)
        endLine = [];
        return;
    end
    endLine = min(endLine, numel(lines));

    % If a tag appears inside the proposed range, cut at the line before first tag
    for i = startLine:endLine
        s = strtrim(lines{i});
        if startsWith(s, '[') && endsWith(s, ']')
            endLine = i - 1;
            break;
        end
    end

    % Remove trailing empty lines
    while endLine >= startLine && isempty(strtrim(lines{endLine}))
        endLine = endLine - 1;
    end

    if endLine < startLine
        endLine = startLine - 1; % empty block
    end
end

function names = local_parseHeader(headerLine)
% Converts tokens like "x [ m ]" to valid MATLAB names like x_m
    raw = strtrim(strsplit(headerLine, ','));
    names = cell(size(raw));

    for i = 1:numel(raw)
        token = strtrim(raw{i});

        m = regexp(token, '^(.*?)\[\s*(.*?)\s*\]\s*$', 'tokens', 'once');
        if ~isempty(m)
            base = strtrim(m{1});
            unit = strtrim(m{2});
        else
            base = token;
            unit = '';
        end

        base = regexprep(base, '[^A-Za-z0-9]+', '_');
        base = regexprep(base, '(^_+|_+$)', '');

        unit = regexprep(unit, '[^A-Za-z0-9]+', '_');
        unit = regexprep(unit, '(^_+|_+$)', '');

        if ~isempty(unit)
            nm = [base '_' unit];
        else
            nm = base;
        end

        if isempty(nm)
            nm = sprintf('Var%d', i);
        end

        names{i} = matlab.lang.makeValidName(nm);
    end
end
