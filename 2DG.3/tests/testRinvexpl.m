classdef testRinvexpl < matlab.unittest.TestCase
    % Tests residual evaluation
    
    methods (Test)
        
        function testConvection(testCase)
            % DRIVER FOR LINEAR CONVECTION PROBLEM                          
            %
            d0=fileparts([pwd,filesep]);
            addpath([d0,'/../app/convection']);   % Add path to application subdirectory

            m      = 2;
            n      = 2;
            porder = 3;

            time  = 2*pi;
            dt    = 2.e-02;
            nstep = 5;
            ncycl = ceil(time/(nstep*dt));

            mesh = mkmesh_square(m,n,porder,1);
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
