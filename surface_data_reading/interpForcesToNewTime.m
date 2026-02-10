function [Fx_new, Fy_new, Fz_new, time_new] = interpForcesToNewTime(time, Fx, Fy, Fz, time_new, method, extrapVal)
%INTERPFORCESTONEWTIME Interpolate force time series (Fx,Fy,Fz) onto time_new.
%
%   [Fx_new, Fy_new, Fz_new, time_new] = interpForcesToNewTime(time, Fx, Fy, Fz, time_new)
%   [Fx_new, Fy_new, Fz_new, time_new] = interpForcesToNewTime(..., method)
%   [Fx_new, Fy_new, Fz_new, time_new] = interpForcesToNewTime(..., method, extrapVal)
%
% Inputs
%   time      : Nx1 or 1xN time vector
%   Fx,Fy,Fz  : force vectors, same length as time (N)
%   time_new  : new time vector to interpolate onto
%   method    : (optional) 'linear' (default), 'pchip', 'spline', 'makima', etc.
%   extrapVal : (optional) extrapolation behavior:
%               - if omitted/empty: no extrapolation; values outside range become NaN
%               - if numeric scalar: use that value outside range
%               - if char/string 'extrap': extrapolate using selected method
%
% Outputs
%   Fx_new, Fy_new, Fz_new : interpolated forces at time_new
%   time_new               : echoed back for convenience

    if nargin < 6 || isempty(method)
        method = 'linear';
    end

    if nargin < 7
        extrapVal = [];
    end

    % Force vectors to column form for consistent indexing
    time = time(:);
    Fx   = Fx(:);
    Fy   = Fy(:);
    Fz   = Fz(:);
    time_new = time_new(:);

    % Basic validation
    n = numel(time);
    if any([numel(Fx), numel(Fy), numel(Fz)] ~= n)
        error('Fx, Fy, and Fz must have the same number of elements as time.');
    end

    % Sort by time to satisfy interp1 requirements
    [time_s, idx] = sort(time);
    Fx_s = Fx(idx);
    Fy_s = Fy(idx);
    Fz_s = Fz(idx);

    % Interpolate
    if isempty(extrapVal)
        Fx_new = interp1(time_s, Fx_s, time_new, method);
        Fy_new = interp1(time_s, Fy_s, time_new, method);
        Fz_new = interp1(time_s, Fz_s, time_new, method);
    else
        Fx_new = interp1(time_s, Fx_s, time_new, method, extrapVal);
        Fy_new = interp1(time_s, Fy_s, time_new, method, extrapVal);
        Fz_new = interp1(time_s, Fz_s, time_new, method, extrapVal);
    end
end
