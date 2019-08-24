function [ mpsOut, norm1, norm2 ] = normalize( mpsIn )
    %NORMALIZE an infinite mps in canonical form as set out in  README
    %   mpsOut = normalize(mpsIn), where mpsIn is an infinite
    %   translationally invariant mps, has normalized Schmidt spectrum and unit 
    %   largest eigenvalue of its transfer matrix. If mpsIn is in canonical 
    %   form, These together imply that it also has a unit length: 
    %   <mpsOut|mpsOut> = 1.
    %
    %   optionally returns the norms, norm1 & norm2, calculated:
    %     lambdaOut = lambdaIn / norm1
    %     T_Out     = T_In     / norm2,
    %   where T=A for left/right canonical form and G for Vidal form.
    
    FORM = mpsIn{1};
    
    % find lamda*lambda':=Lambda and relevant transfer matrix
    
    switch FORM
        case {1,2}
            Lambda = mpsIn{3};
            A = mpsIn{2};
            E = ncon({A,conj(A)}, {[-1 -2 1], [-3 -4 1]});
        case 3
            lambda = mpsIn{3};
            Lambda = lambda*lambda';
            G = mpsIn{2}; % \Gamma tensor
            E = ncon({G, lambda, conj(G), conj(lambda)},... % this is usually
                     {[-1 1 2], [1 -2], [-3 3 2], [3 -4]}); % called 'R'
        otherwise
            err.identifier = 'normalize:InvalidFormInput';
            err.message = 'Invalid value of FORM = mpsIn{1}';
            error(err);
    end
    
    % normalize and output as a cell array
    % Condition I. tr(lambda.lambda') = 1
    norm1 = sqrt(trace(Lambda));
    
    % Condition II. unit eigenvalues for canonical form conditions
    IEI = ncon({E}, {[1 2 1 2]}); % with E=R if we're in Vidal form
    D = size(E,1);
    norm2 = sqrt(IEI/D);
    
    switch FORM
        case {1,2}
            mpsOut = {FORM, A/norm2, Lambda/(norm1^2)};
        case 3
            mpsOut = {FORM, G/norm2, lambda/norm1};
    end
end