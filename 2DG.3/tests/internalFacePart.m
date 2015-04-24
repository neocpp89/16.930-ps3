function [r, myr] = internalFacePart(m, n, porder)
    mesh = mkmesh_square(m,n,porder,1);
    mesh = mkmesh_distort(mesh);
    master = mkmaster(mesh,2*porder);
    app = mkapp;
    app.finvv = @zv;
    app.finvb = @(up,np,ib,ui,p,param,time) zeros(size(up));

    app.bcm = [1,1,1,1];                    % Manually set boundary conditions
    app.bcs = [0];

    vf = @(p) fliplr(p-0.5)*diag([-1,1]);   % Rotating field
    app.arg = {vf};

    init  = @(dg) exp(-120*((dg(:,1,:)-0.6).^2 + (dg(:,2,:)-0.5).^2));  % Gaussian hill
    u = initu(mesh,app,{init});

    tm = 0;
    r = rinvexpl(master,mesh,app,u,tm);
    disp('-----');
    myr = myrinvexpl(master,mesh,app,u,tm);
end
