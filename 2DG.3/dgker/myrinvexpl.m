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

npl = size(mesh.plocal, 1);
nc = app.nc;

% loop through elements
nt = size(mesh.t, 1);
for i=1:nt
    p(:,:) = mesh.dgnodes(:,:,i);
    uv = u(:,:,i);
    r(:,:,i) = r(:,:,i) + app.finvv(uv, p, app.arg, time);
end

% loop through faces
nf = size(mesh.f, 1);
for i=1:nf
end

ni = sum(mesh.f(:,4) > 0);
% loop through internal faces


% divide by mass matrix
for k=1:size(mesh.t,1)
    for j=1:app.nc
        r(:,j,k) = master.mass \ r(:,j,k);
    end
end
end