%DRIVER FOR THE EULER EQUATIONS IN A CHANNEL (WITH OBSTACLE)
%
d0=fileparts([pwd,filesep]);
addpath([d0,'/Euler']);       % Add path to application subdirectory

p = [1 3];
xm = [0 21];
xn = [0 12];

xm(1) = ceil(sqrt(3)*xm(2));
xn(1) = ceil(sqrt(3)*xn(2));

meshes = cell(numel(p), 1);
solutions = cell(numel(p), 1);

parfor j=1:numel(p)
    m      = xm(j);
    n      = xn(j);
    porder = p(j);

    time  = 10;
    dt    = 2.e-03;
    nstep = 5;
    ncycl = ceil(time/(nstep*dt));
    gam = 1.4;
    Minf = 0.3;                    % Infinity conditions
    ui = [ 1, 1, 0, 0.5+1/(gam*(gam-1)*Minf^2)];

    db   = 0.2;                     % Height of channel constriction
    mesh = mkmesh_square(m,n,porder,1);
    mesh = mkmesh_duct(mesh, db,0.0,1);
    master = mkmaster(mesh,2*porder);
    app = mkapp;

    app.bcm = [2,1,2,1];  % 2: Wall, 1: Far-field, Bottom-Top: Wall, Left-Right: Far-field
    app.bcs = [ui; 0*ui];
    app.arg = {gam};

    vv = @(dg) ones(size(dg,1),size(dg,3))*ui(3);

    rho  = @(dg) 1+0.01*exp(-80*((dg(:,1,:)-1.5).^2 + (dg(:,2,:)-0.5).^2));
    ru   = @(dg) rho(dg)*ui(2)/ui(1);
    rv   = @(dg) rho(dg)*ui(3)/ui(1);
    pinf = 1/(gam*Minf^2);
    rE   = @(dg) 0.5*(ru(dg).^2+rv(dg).^2)./rho(dg) + pinf/(gam-1);

    u = initu(mesh,app,{rho,ru,rv,rE});
    tm = 0.0;
    
    tic; rk4(@myrinvexpl,master,mesh,app,u,tm,dt,nstep); s = toc;
    fprintf('This run should take about %d seconds.\n', s*ncycl);

    for i = 1:ncycl    
        % scaplot(mesh,eulereval(u,'r',gam),[],[],1); axis off;
        
        u = rk4(@myrinvexpl,master,mesh,app,u,tm,dt,nstep);
        tm = tm + nstep*dt;
        fprintf('[%d] step = %d/%d\n', j, i, ncycl);
    end
    meshes{j} = mesh;
    solutions{j} = u;
end
save('channel_obstacle.mat');
