function out = convertForcesnTorqueToLocalCoor(out,R)

Np = numel(out.PatchIDs);
Nt = size(out.F,1);

out.F_local = zeros(Nt,3,Np);
out.T_local = zeros(Nt,3,Np);
out.T_aboutZ = zeros(Nt,Np);   % scalar torque about local z-axis

for k = 1:Np
    Fg = out.F(:,:,k);  % Nt x 3
    Tg = out.T(:,:,k);  % Nt x 3

    % Convert: v_local = (R' * v_global')'
    out.F_local(:,:,k) = (R' * Fg.').';
    out.T_local(:,:,k) = (R' * Tg.').';
    
    % Torque about the defined axis (local z component)
    out.T_aboutX(:,k) = out.T_local(:,1,k);
    out.T_aboutY(:,k) = out.T_local(:,2,k);
    out.T_aboutZ(:,k) = out.T_local(:,3,k);

    out.F_aboutX(:,k) = out.F_local(:,1,k);
    out.F_aboutY(:,k) = out.F_local(:,2,k);
    out.F_aboutZ(:,k) = out.F_local(:,3,k);    
end

end