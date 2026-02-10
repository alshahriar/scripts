function [R, O, ex, ey, ez] = coordFrameFromAxisPoints(O, Zp, XZp)
%COORDFRAMEFROMAXISPOINTS Build right-handed orthonormal frame from axis points.
%
% Inputs: O, Zp, XZp are 1x3 vectors in global coordinates.
% Output:
%   ex,ey,ez : unit basis vectors in GLOBAL coordinates
%   R        : 3x3 rotation matrix whose columns are [ex ey ez]
%             (maps local->global: v_global = R * v_local)

O  = O(:).';  Zp = Zp(:).';  XZp = XZp(:).';

% z-axis
ez = Zp - O;
nz = norm(ez);
if nz < eps, error('Origin and Z Axis Point are identical (cannot define z-axis).'); end
ez = ez / nz;

% x-axis: take vector to XZ-plane point, remove its component along ez
v  = XZp - O;
vx = v - dot(v, ez)*ez;
nx = norm(vx);
if nx < eps
    error('XZ Plane Pt is collinear with z-axis (cannot define x-axis). Choose a different XZ Plane Pt.');
end
ex = vx / nx;

% y-axis (right-handed)
ey = cross(ez, ex);
ny = norm(ey);
if ny < eps, error('Degenerate frame (check points).'); end
ey = ey / ny;

% Re-orthogonalize x to be safe
ex = cross(ey, ez);

R = [ex(:), ey(:), ez(:)];
end
