classdef testRinvexpl < matlab.unittest.TestCase
    % Tests residual evaluation
    
    methods (Test)
        function testVolumePart(testCase)
            m      = 2;
            n      = 2;
            porder = 3;

            mesh = mkmesh_square(m,n,porder,1);
            mesh = mkmesh_distort(mesh);
            master = mkmaster(mesh,2*porder);
            app = mkapp;
            app.finvi = @(up,um,np,p,param,time) zeros(size(up));
            app.finvb = @(up,np,ib,ui,p,param,time) zeros(size(up));

            app.bcm = [1,1,1,1];                    % Manually set boundary conditions
            app.bcs = [0];

            vf = @(p) fliplr(p-0.5)*diag([-1,1]);   % Rotating field
            app.arg = {vf};

            init  = @(dg) exp(-120*((dg(:,1,:)-0.6).^2 + (dg(:,2,:)-0.5).^2));  % Gaussian hill
            u = initu(mesh,app,{init});

            tm = 0;
            r_expected = rinvexpl(master,mesh,app,u,tm);
            disp('-----');
            r_actual = myrinvexpl(master,mesh,app,u,tm);
            testCase.verifyEqual(r_actual, r_expected, 'Abstol', 1e-10);
        end
        function testFull(testCase)
            m      = 2;
            n      = 2;
            porder = 3;

            mesh = mkmesh_square(m,n,porder,1);
            mesh = mkmesh_distort(mesh);
            master = mkmaster(mesh,2*porder);
            app = mkapp;

            app.bcm = [1,1,1,1];                    % Manually set boundary conditions
            app.bcs = [0];

            vf = @(p) fliplr(p-0.5)*diag([-1,1]);   % Rotating field
            app.arg = {vf};

            init  = @(dg) exp(-120*((dg(:,1,:)-0.6).^2 + (dg(:,2,:)-0.5).^2));  % Gaussian hill
            u = initu(mesh,app,{init});

            tm = 0;
            r_expected = rinvexpl(master,mesh,app,u,tm);
            disp('-----');
            r_actual = myrinvexpl(master,mesh,app,u,tm);
            testCase.verifyEqual(r_actual, r_expected, 'Abstol', 1e-10);
        end
    end
end
