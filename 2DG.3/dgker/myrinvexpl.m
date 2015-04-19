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

for i=1:size(mesh.dgnodes, 1)
    for t=1:size(mesh.dgnodes, 3)
        J = [xxi(i,t), xet(i,t); yxi(i,t), yet(i,t)];
        jacobian(:,:,i,t) = J;
        detJ(i,t) = det(J);
    end
end

phi1d(:,:) = master.sh1d(:,1,:);
dphi1d(:,:) = master.sh1d(:,2,:);
phi(:,:) = master.shap(:,1,:);
dphidxi(:,:) = master.shap(:,2,:);
dphideta(:,:) = master.shap(:,3,:);

npl = size(mesh.plocal, 1);
nc = app.nc;

% loop through all faces
nf = size(mesh.f, 1);
for i=1:nf
end

% loop through internal faces
ni = sum(mesh.f(:,4) > 0);
for i=1:ni
end

% loop through elements
nt = size(mesh.t, 1);
for i=1:nt
    p(:,:) = mesh.dgnodes(:,:,i);
    uv(:,:) = u(:,:,i);
    pg = phi'*p;
    ug = phi'*uv;
    finvv = app.finvv(ug, pg, app.arg, time);
    dphidx = dphidxi*pg
    dphidy = dphideta*pg
    size(finvv)
    size(jcw(:,i))
    r(:,:,i) = r(:,:,i) + jcw(:, i).*finvv;
end


% divide by mass matrix
for k=1:size(mesh.t,1)
    for j=1:app.nc
        r(:,j,k) = master.mass \ r(:,j,k);
    end
end
end
