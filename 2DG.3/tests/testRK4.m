classdef testRK4 < matlab.unittest.TestCase
    % Tests residual evaluation with rk4

    methods (Test)
        function testRK4Full(testCase)
            m = 2;
            n = 2;
            porder = 3;
            
            mesh = mkmesh_square(m,n,porder,1);
            mesh = mkmesh_distort(mesh);
            master = mkmaster(mesh,2*porder);
            app = mkapp;

            app.bcm = [1,1,1,1];                    % Manually set boundary conditions
            app.bcs = [0];

            vf = @(p) fliplr(p-0.5)*diag([-1,1]);   % Rotating field
            app.arg = {vf};
            app.arg = { @(p) ones(size(p)) };

            init  = @(dg) exp(-120*((dg(:,1,:)-0.6).^2 + (dg(:,2,:)-0.5).^2));  % Gaussian hill
            u = initu(mesh,app,{init});

            tm = 0;
            time = tm;
            dt = 2e-2;
            nstep = 5;
            
            k1mr = dt*myrinvexpl( master, mesh, app, u       , time       );
            k1r = dt*rinvexpl( master, mesh, app, u       , time       );
            
            k2r = dt*rinvexpl( master, mesh, app, u+0.5*k1r, time+0.5*dt);
            k2mr = dt*myrinvexpl( master, mesh, app, u+0.5*k1r, time+0.5*dt);
            
            k3r = dt*rinvexpl( master, mesh, app, u+0.5*k2r, time+0.5*dt);
            k3mr = dt*myrinvexpl( master, mesh, app, u+0.5*k2mr, time+0.5*dt);
            
            k4r = dt*rinvexpl( master, mesh, app, u+    k3r, time+    dt);
            k4mr = dt*myrinvexpl( master, mesh, app, u+    k3mr, time+    dt);

            testCase.verifyEqual(k1mr, k1r, 'Abstol', 1e-10);
            testCase.verifyEqual(k2mr, k2r, 'Abstol', 1e-10);
            testCase.verifyEqual(k3mr, k3r, 'Abstol', 1e-10);
            testCase.verifyEqual(k4mr, k4r, 'Abstol', 1e-10);

            u_actual = rk4(@myrinvexpl,master,mesh,app,u,tm,dt,nstep);
            u_expected = rk4(@rinvexpl,master,mesh,app,u,tm,dt,nstep);

            testCase.verifyEqual(u_actual, u_expected, 'Abstol', 1e-10);
        end
    end
end
