function r = myrinvexpl(master,mesh,app,u,time)
%RINVEXPL Calculates the residual vector for explicit time stepping.
%   R=RINVEXPL(MASTER,MESH,APP,U,TIME)
%
%      MASTER:       Master structure
%      MESH:         Mesh structure
%      APP:          Application structure
%      U(NPL,NC,NT): Vector of unknowns
%                    NPL = size(mesh.plocal,1)
%                    NC = app.nc (number of equations in system)
%                    NT = size(mesh.t,1)
%      TIME:         Time
%      R(NPL,NC,NT): Residual vector (=dU/dt) (already divided by mass
%                    matrix)                             
r = zeros(size(u));

xxi = squeeze(master.shap(:,2,:))'*squeeze(mesh.dgnodes(:,1,:));
xet = squeeze(master.shap(:,3,:))'*squeeze(mesh.dgnodes(:,1,:));
yxi = squeeze(master.shap(:,2,:))'*squeeze(mesh.dgnodes(:,2,:));
yet = squeeze(master.shap(:,3,:))'*squeeze(mesh.dgnodes(:,2,:));
detJ = (xxi.*yet - xet.*yxi);

phi1d(:,:) = master.sh1d(:,1,:);
dphi1d(:,:) = master.sh1d(:,2,:);
phi(:,:) = master.shap(:,1,:);
dphidxi(:,:) = master.shap(:,2,:);
dphideta(:,:) = master.shap(:,3,:);

% loop through all faces
nf = size(mesh.f, 1);
for i=1:nf
    kl = mesh.f(i,3); % left element index
    eli = find(abs(mesh.t2f(kl,:)) == i); % left element edge index
    elo = 1 + (mesh.t2f(kl,eli) < 0); % orientation of edge on left element

    ep = mesh.dgnodes(master.perm(:, eli, elo), :, kl);
    epg = phi1d'*ep; % location of gauss points on edge

    % tangent vector
    xsg = dphi1d'*ep(:,1);
    ysg = dphi1d'*ep(:,2);
    dsg = sqrt(ysg.*ysg + xsg.*xsg);

    % outward normal
    nepg = [ysg ./ dsg, -xsg ./ dsg];

    ul = u(master.perm(:, eli, elo), :, kl); % solution on left element's edge
    ulg = phi1d'*ul; % and at gauss points on edge

    scale = (master.gw1d .* sqrt(xsg.*xsg + ysg.*ysg));
    S = diag(scale);

    kr = mesh.f(i,4); % right element (tells if internal face or boundary)
    if (kr >= 0)
        % internal face, get right element quantities
        eri = find(abs(mesh.t2f(kr,:)) == i);
        ero = 1 + (mesh.t2f(kr,eri) < 0);
        ur = u(master.perm(:, eri, ero), :, kr);
        urg = phi1d'*ur;

        lfinvi = app.finvi( ulg, urg, nepg, epg, app.arg, time);
        rfinvi = app.finvi( urg, ulg, -nepg, epg, app.arg, time);

        r(master.perm(:, eli, elo), :, kl) = r(master.perm(:, eli, elo), :, kl) - (phi1d * S * lfinvi);
        r(master.perm(:, eri, ero), :, kr) = r(master.perm(:, eri, ero), :, kr) - (phi1d * S * rfinvi);
    else
        % external face
        ibt = app.bcm(-kr);
        finvb = app.finvb( ulg, nepg, ibt, app.bcs(ibt, :), epg, app.arg, time);
        r(master.perm(:, eli, elo), :, kl) = r(master.perm(:, eli, elo), :, kl) - (phi1d * S * finvb);
    end
end


% loop through elements
nt = size(mesh.t, 1);
for i=1:nt
    p(:,:) = mesh.dgnodes(:,:,i);
    uv(:,:) = u(:,:,i);
    pg = phi'*p;
    ug = phi'*uv;

    % components of inverse jacobian
    xix = diag(yet(:,i) ./ detJ(:,i));
    etax = diag(-yxi(:,i) ./ detJ(:,i));
    xiy = diag(-xet(:,i) ./ detJ(:,i));
    etay = diag(xxi(:,i) ./ detJ(:,i));
    scale = (master.gwgh .* detJ(:,i));
    S = diag(scale);

    % derivatives with respect to global coordinates
    dphidx = (dphidxi*xix + dphideta*etax);
    dphidy = (dphidxi*xiy + dphideta*etay);

    % mass matrix if non-constant jacobians
    % mm = phi*diag(master.gwgh .* detJ(:, i))*phi';

    [fvx, fvy] = app.finvv(ug, pg, app.arg, time);
    finvv = dphidx * S * fvx + dphidy * S * fvy;
    r(:,:,i) = r(:,:,i) + finvv;
    % r(:,:,i) = r(:,:,i) + (mm \ finvv);
end

% divide by mass matrix
nt = size(mesh.t, 1);
for i=1:nt
    % mass matrix if non-constant jacobians
    mm = phi*diag(master.gwgh .* detJ(:, i))*phi';
    r(:,:,i) = (mm \ r(:,:,i));
end

end
