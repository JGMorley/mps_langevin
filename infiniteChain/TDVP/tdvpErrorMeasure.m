function epsilon = tdvpErrorMeasure(A, H)
    %% Calculate epsilon, the norm of the part of (H-<H>)|psi> that is lost
    %  on projection to the tangent space 
    %  As defined in Haegeman et al, PRB 88 075133 (2013)
    %  DOI:PhysRevB.88.075133
    
    % 1. find sqrtl, sqrtr, VL, VR
    [ ~, right, left] = normalizeMPS(A);
    sqrtl = squareDecomposition(left);
    sqrtr = squareDecomposition(right);
    
    [~,D,d] = size(A);
    VL = leftGaugeTangentBasis( sqrtl, A, d, D);
    VR = rightGaugeTangentBasis(sqrtr, A, d, D);
    
    % 2. compute Y and epsilon   
    Y = calculateY(A,sqrtl,sqrtr,VL,VR,H);
    epsilon = sqrt(trace(Y*Y'));
end