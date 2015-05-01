function [logical_idx, polar] = find_nodes_at_radius(m, r, tol)
% m is the mesh
% r is the radius you want
% tol is the tolerance in radius

% returns:
%   logical_idx - logical index of dgnodes which have required radius
%   polar - dgnode coordinates converted to polar
    polar = topolar(m.dgnodes);
    rad(:,:) = polar(:,1,:);
    logical_idx = (abs(rad - r) < tol);
end

function [polarnodes] = topolar(dgnodes)
    polarnodes(:,1,:) = sqrt(dgnodes(:,1,:).^2 + dgnodes(:,2,:).^2);
    polarnodes(:,2,:) = atan2(dgnodes(:,2,:), dgnodes(:,1,:));
end
