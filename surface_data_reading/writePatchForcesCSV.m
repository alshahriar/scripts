function T = writePatchForcesCSV(csvFile, t, out, patchNames, varargin)
%WRITEPATCHFORCESCSV Write per-patch forces (Fx,Fy,Fz) vs time to a CSV file.
%
% T = writePatchForcesCSV(csvFile, t, out, patchNames, 'Name',Value,...)
%
% Inputs
%   csvFile    : output CSV filename (e.g., 'patch_forces.csv')
%   t          : Nt x 1 time vector
%   out        : struct from pressureForcesByPatches(), must contain:
%                out.PatchIDs (Np x 1), out.Fx/out.Fy/out.Fz (Nt x Np)
%   patchNames : Np x 1 cellstr/string array of names for each patch column header,
%                in the SAME order as out.PatchIDs. Example: ["inlet","shroud","hub"]
%
% Name-value options
%   'IncludePatchIDsRow' : true/false (default false)  % writes a 2nd header row with PatchIDs
%   'Delimiter'          : ',' (default) or other for writetable
%
% Output
%   T : table written to disk

p = inputParser;
p.addRequired('csvFile', @(s) ischar(s) || isstring(s));
p.addRequired('t', @(x) isnumeric(x) && isvector(x));
p.addRequired('out', @(s) isstruct(s) && isfield(s,'PatchIDs') && isfield(s,'Fx') && isfield(s,'Fy') && isfield(s,'Fz'));
p.addRequired('patchNames', @(x) iscellstr(x) || isstring(x));
p.addParameter('IncludePatchIDsRow', false, @(x) islogical(x) && isscalar(x));
p.addParameter('Delimiter', ',', @(x) ischar(x) || isstring(x));
p.parse(csvFile, t, out, patchNames, varargin{:});

t = t(:);
Nt = numel(t);

PatchIDs = out.PatchIDs(:);
Np = numel(PatchIDs);

% Normalize patchNames
patchNames = string(patchNames(:));
if numel(patchNames) ~= Np
    error('patchNames must have %d entries (same as number of patches in out.PatchIDs).', Np);
end

% Validate sizes
if size(out.Fx,1) ~= Nt || size(out.Fy,1) ~= Nt || size(out.Fz,1) ~= Nt
    error('t length (%d) must match number of time steps in out.Fx/out.Fy/out.Fz.', Nt);
end
if size(out.Fx,2) ~= Np || size(out.Fy,2) ~= Np || size(out.Fz,2) ~= Np
    error('out.Fx/out.Fy/out.Fz must have %d columns (one per patch).', Np);
end

% Make valid, unique variable names for table columns
baseNames = matlab.lang.makeValidName(patchNames);
baseNames = matlab.lang.makeUniqueStrings(baseNames);

% Build table
T = table(t, 'VariableNames', {'t'});

for k = 1:Np
    T.("Fx_" + baseNames(k)) = out.F_aboutX(:,k);
end
for k = 1:Np
    T.("Fy_" + baseNames(k)) = out.F_aboutY(:,k);
end
for k = 1:Np
    T.("Fz_" + baseNames(k)) = out.F_aboutZ(:,k);
end

for k = 1:Np
    T.("Tx_" + baseNames(k)) = out.T_aboutX(:,k);
end
for k = 1:Np
    T.("Ty_" + baseNames(k)) = out.T_aboutY(:,k);
end
for k = 1:Np
    T.("Tz_" + baseNames(k)) = out.T_aboutZ(:,k);
end


% Write
writetable(T, csvFile, 'Delimiter', p.Results.Delimiter);

% Optional: add a second header row with PatchIDs (Excel-friendly trick)
if p.Results.IncludePatchIDsRow
    % Read the file back as text and prepend a row. This keeps writetable simplicity.
    header1 = string(T.Properties.VariableNames);
    header2 = strings(size(header1));
    header2(1) = "PatchID";
    % Map IDs into the Fx/Fy/Fz columns by patch
    % Columns are: t | Fx_* (Np) | Fy_* (Np) | Fz_* (Np)
    header2(2:1+Np)         = string(PatchIDs);
    header2(2+Np:1+2*Np)    = string(PatchIDs);
    header2(2+2*Np:1+3*Np)  = string(PatchIDs);

    txt = fileread(csvFile);
    newTxt = strjoin(header1, p.Results.Delimiter) + newline + ...
             strjoin(header2, p.Results.Delimiter) + newline + ...
             txt; % txt already includes header+data; we will remove first header line
    % Remove the original header line from txt:
    firstNL = find(txt==newline, 1, 'first');
    if ~isempty(firstNL)
        txtNoHeader = txt(firstNL+1:end);
        newTxt = strjoin(header1, p.Results.Delimiter) + newline + ...
                 strjoin(header2, p.Results.Delimiter) + newline + ...
                 txtNoHeader;
    end

    fid = fopen(csvFile, 'w');
    fwrite(fid, newTxt);
    fclose(fid);
end
end
