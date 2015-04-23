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

for i=1:size(xxi, 1)
    for t=1:size(mesh.dgnodes, 3)
        J = [xxi(i,t), xet(i,t); yxi(i,t), yet(i,t)];
        jacobian(:,:,i,t) = J;
        detJ = 1.0 ./ (xxi.*yet - xet.*yxi);
        jacobian_inverse(:,:,i,t) = J^(-1);
    end
end

phi1d(:,:) = master.sh1d(:,1,:);
dphi1d(:,:) = master.sh1d(:,2,:);
phi(:,:) = master.shap(:,1,:);
dphidxi(:,:) = master.shap(:,2,:);
dphideta(:,:) = master.shap(:,3,:);

dsdxi = sqrt(yxi.^2 + xxi.^2);

npl = size(mesh.plocal, 1);
nc = app.nc;

% loop through faces of elements
nt = size(mesh.t, 1);
for it=1:nt
    for ie=1:size(mesh.t2f, 2)
        orientation = 1+(mesh.t2f(it, ie) < 0);
        ep = squeeze(mesh.dgnodes(master.perm(:,ie,orientation)));
        epg = phi1d'*ep; 
    end
end

% loop through internal faces
ni = sum(mesh.f(:,4) > 0);
for i=1:ni
    epl = mesh.f(i,1) % left point index 
    epr = mesh.f(i,2)
    kl = mesh.f(i,3) % left element index
    kr = mesh.f(i,4)
    eli = find(abs(mesh.t2f(kl,:)) == i) % left element edge index
    eri = find(abs(mesh.t2f(kr,:)) == i)
    elo = 1 + (mesh.t2f(kl,eli) < 0) % orientation of edge on left element
    ero = 1 + (mesh.t2f(kr,eri) < 0)
end

% loop through elements
nt = size(mesh.t, 1);
for i=1:nt
    p(:,:) = mesh.dgnodes(:,:,i);
    uv(:,:) = u(:,:,i);
    pg = phi'*p;
    ug = phi'*uv;
    xix = diag(yet(:,i) ./ detJ(:,i));
    etax = diag(-yxi(:,i) ./ detJ(:,i));
    xiy = diag(-xet(:,i) ./ detJ(:,i));
    etay = diag(xxi(:,i) ./ detJ(:,i));
    scale = diag(master.gwgh .* detJ(:,i));
    dphidx = (dphidxi*xix + dphideta*etax);
    dphidy = (dphidxi*xiy + dphideta*etay);

    finvv = app.finvv(ug, pg, app.arg, time);
    r(:,:,i) = r(:,:,i) + finvv;
end


% divide by mass matrix
for k=1:size(mesh.t,1)
    for j=1:app.nc
        r(:,j,k) = master.mass \ r(:,j,k);
    end
end
end
