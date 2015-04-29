%DRIVER FOR LINEAR CONVECTION PROBLEM                          
%
d0=fileparts([pwd,filesep]);
addpath([d0,'/convection']);   % Add path to application subdirectory

n = [4 6 8 12 16];
p = [2 3 4];

figure;
hold all;
for pp=p
    err = [];
    for x=n
        m      = x;                    
        n      = x;
        porder = pp;

        time  = 2*pi;
        dt    = 1.e-02;
        nstep = 5;
        ncycl = ceil(time/(nstep*dt));

        mesh = mkmesh_square(m,n,porder,1);
        % mesh = mkmesh_distort(mesh,0.05);     % Uncomment for mesh distortion 
        master = mkmaster(mesh,2*porder);
        app = mkapp;

        app.bcm = [1,1,1,1];                    % Manually set boundary conditions
        app.bcs = [0];

        vf = @(p) fliplr(p-0.5)*diag([-1,1]);   % Rotating field
        app.arg = {vf};

        init  = @(dg) exp(-120*((dg(:,1,:)-0.6).^2 + (dg(:,2,:)-0.5).^2));  % Gaussian hill
        u = initu(mesh,app,{init});
        u0 = u;

        tm = 0.0;
        for i=1:ncycl   
            % scaplot(mesh,u,[-0.1,1.2],[],1); axis off; 
            u = rk4(@myrinvexpl,master,mesh,app,u,tm,dt,nstep);
            tm = tm + nstep*dt;
        end

        while tm < time                         % Complete an exact revolution
            u = rk4(@myrinvexpl,master,mesh,app,u,tm,dt,1);
            tm = tm + dt;
        end
        % scaplot(mesh,u,[-0.1,1.2],[],1); axis off;

        xxi = squeeze(master.shap(:,2,:))'*squeeze(mesh.dgnodes(:,1,:));
        xet = squeeze(master.shap(:,3,:))'*squeeze(mesh.dgnodes(:,1,:));
        yxi = squeeze(master.shap(:,2,:))'*squeeze(mesh.dgnodes(:,2,:));
        yet = squeeze(master.shap(:,3,:))'*squeeze(mesh.dgnodes(:,2,:));
        detJ = (xxi.*yet - xet.*yxi);

        phi = squeeze(master.shap(:,1,:));

        nt = size(mesh.t, 1);
        errsq = 0;
        for i=1:nt
            uv = squeeze((u(:,:,i) - u0(:,:,i)) .* (u(:,:,i) - u0(:,:,i)));
            ug = phi'*uv;

            scale = (master.gwgh .* detJ(:,i));
            S = diag(scale);

            errsq = errsq + sum(S * ug);
        end
        err = [err sqrt(errsq)];
        disp(x);
    end
    plot(1./n, err);
    legend(sprintf('Order %d', pp));
end
hold off;

