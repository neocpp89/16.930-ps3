classdef testRinvexpl < matlab.unittest.TestCase
    % Tests residual evaluation

    methods (Test)
        function testSmallDomainVolumePart(testCase)
            [r_expected, r_actual] = volumePart(2, 2, 3);
            testCase.verifyEqual(r_actual, r_expected, 'Abstol', 1e-10);
        end
        function testSmallDomainInternalFacePart(testCase)
            [r_expected, r_actual] = internalFacePart(2, 2, 3);
            testCase.verifyEqual(r_actual, r_expected, 'Abstol', 1e-10);
        end
        function testSmallDomainBoundaryFacePart(testCase)
            [r_expected, r_actual] = boundaryFacePart(2, 2, 3);
            testCase.verifyEqual(r_actual, r_expected, 'Abstol', 1e-10);
        end
        function testSmallDomainFullResidual(testCase)
            [r_expected, r_actual] = fullResidual(2, 2, 3);
            testCase.verifyEqual(r_actual, r_expected, 'Abstol', 1e-10);
        end
        function testVolumePart(testCase)
            for p=1:3
                for m=2:2:6
                    for n=2:2:6
                        [r_expected, r_actual] = volumePart(m, n, p);
                        testCase.verifyEqual(r_actual, r_expected, 'Abstol', 1e-10);
                    end
                end
            end
        end
        function testInternalFacePart(testCase)
            for p=1:3
                for m=2:2:6
                    for n=2:2:6
                        [r_expected, r_actual] = internalFacePart(m, n, p);
                        testCase.verifyEqual(r_actual, r_expected, 'Abstol', 1e-10);
                    end
                end
            end
        end
        function testBoundaryFacePart(testCase)
            for p=1:3
                for m=2:2:6
                    for n=2:2:6
                        [r_expected, r_actual] = boundaryFacePart(m, n, p);
                        testCase.verifyEqual(r_actual, r_expected, 'Abstol', 1e-10);
                    end
                end
            end
        end
        function testFull(testCase)
            for p=1:3
                for m=2:2:6
                    for n=2:2:6
                        [r_expected, r_actual] = fullResidual(m, n, p);
                        testCase.verifyEqual(r_actual, r_expected, 'Abstol', 1e-10);
                    end
                end
            end
        end
    end
end
