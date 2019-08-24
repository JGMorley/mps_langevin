function mpsOut = mpsGenerator(p, D, FORM)
    %Returns an MPS with physical index and bond index dimensions equal to p, D
    %  mpsOut = mspGenerator(p, D, FORM) returns a cell array mpsOut where
    %    for FORM = 3,      mpsOut is a randomized MPS in Vidal form,
    %    for FORM = 1 or 2, mpsOut is a randomized MPS in terms of 'A' tensors
    %
    %  see README.md for more on MPS storage conventions
    
    mpsOut{1} = FORM;
    GorA = rand([D,D,p]) + rand([D,D,p])*1i;
    
    lambda = diag(rand([D,1])); % Does this need to be diagonal?
    lambda = rand([D,D]) + rand([D,D])*1i;
    
    switch FORM
        case 1 % A = lambda.G
            mpsOut{2} = ncon({lambda,GorA},{[-1 1], [1 -2 -3]});
            mpsOut{3} = ncon({lambda, conj(lambda)}, {[1 -1], [1 -2]});
        case 2 % A = G.lambda
            mpsOut{2} = ncon({GorA,lambda},{[-1 1 -3], [1 -2]});
            mpsOut{3} = ncon({lambda, conj(lambda)}, {[-1 1], [-2 1]});
        case 3 % Vidal form
            mpsOut{2} = GorA;
            mpsOut{3} = lambda;
    end
end