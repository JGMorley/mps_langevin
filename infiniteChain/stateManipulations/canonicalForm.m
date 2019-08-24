function [mpsOut, A, normFactor] = canonicalForm(mpsIn, TOL)
    %Puts mpsIn into canonical form, using tolerance of TOL digits for checks
    %  mpsOut = canonicalForm(mpsIn, TOL) returns an MPS specified by a
    %    cell-array mpsOut such that each bond index corresponds to a
    %    Schmidt decomposition
    %  NB also normalizes the state so that 
    %     I. tr(lambda.lambda') = 1, and
    %     II. eigenvalues for canonical form conditions are unity.
    %
    %  optional output normFactor gives the overall factor that the state has 
    %  been effectively multiplied by in order to normalize  
    %
    %  see eg. DOI: 10.1103/PhysRevB.78.155117
    
    if nargin < 2
        TOL = 10;
    end
    
    FORM = mpsIn{1};
    
    switch FORM
        case 1 % left canonical form
            [mpsOut, A, normFactor] = AToCanonicalForm(mpsIn, TOL);
        case 2 % right canonical form
            [mpsOut, A, normFactor] = AToCanonicalForm(mpsIn, TOL); 
        case 3 % Vidal form - A is Gamma here
            [mpsOut, A, normFactor] = GLToCanonicalForm(mpsIn, TOL); 
    end
end

function [mpsOut, gammaOut, normFactor] = GLToCanonicalForm(mpsIn, TOL)
    % Transforms an mpsIn given in Vidal form to Vidal canonical form
    
    if nargin < 2
        TOL = 10;
    end
    
    FORM     = mpsIn{1};
    gammaIn  = mpsIn{2};
    lambdaIn = mpsIn{3};
    
    if FORM ~= 3
        error('called function with mpsIn of wrong type (FORM != 0)')
    end  
  
    % extract degree of bond index (D)
    D = size(gammaIn,1);
    
    % Find R and L matrices as matrices with left and right indices grouped
    L = ncon({gammaIn, lambdaIn, conj(gammaIn), conj(lambdaIn)},...
             {[2 -2 1], [-1 2], [3 -4 1], [-3 3]});
    L = permute(L,[3,1,4,2]);
    L = reshape(L,[D^2,D^2]);
    
    R = ncon({gammaIn, lambdaIn, conj(gammaIn), conj(lambdaIn)},...
             {[-1 2 1], [2 -2], [-3 3 1], [3 -4]});
    R = permute(R,[3,1,4,2]);
    R = reshape(R,[D^2,D^2]);
 
    % Find dominant eigenvectors and eigenvalues of R and L
    [rightV, righteta] = eigs(R,1);   
    [leftV,  lefteta]  = eigs(L.',1); % NB transpose used for left eigenvector
    
    if lefteta - righteta > 1e-10
        error('Left and right dominant eigenvalues unequal');
    end
    if not(round(imag(righteta),TOL) == 0)
        error('eigenvalue not real within tolerance')
    end
    
    % Decompose leftV, rightV into squares to find X, Y by reshaping,
    % correcting sign and using squareDecomposition()
    rightV = reshape(rightV, [D, D]).'; % .' to account for reshape()
    leftV = reshape(leftV, [D, D]);     % transposed following Orus2009 fig 2
   
    % normalize by dividing by phase of first element
    rightV = rightV/sign(rightV(1,1)); % sign(z) = z/abs(z) for complex z
    leftV = leftV/sign(leftV(1,1));
    
    % [rightV, leftV] should be Hermitian and non-negative     
    X = squareDecomposition(rightV,TOL); % rightV = X*X'
    Y = squareDecomposition(leftV, TOL)';% leftV  = Y'*Y 
                                         % (leftV.'*L = eta leftV.')
    
    % Perform SVD on Y.lambdaIn.X to find U.lambdaOut.V'
    [U, lambdaOut, V] = svd(Y * lambdaIn * X);
    V = V'; % so that U.lambdaOut.V = Y.lambdaIn.X
    
    % Calculate gammaOut
    gammaOut = ncon({V, inv(X), gammaIn, inv(Y), U},...
                    {[-1 1], [1 2], [2 3 -3], [3 4], [4 -2]});
    
    % normalize and output as a cell array
    % Condition I. tr(lambda.lambda') = 1
    norm1 = sqrt(trace(lambdaOut*lambdaOut'));
    lambdaOut = lambdaOut/norm1;
    
    % Condition II. unit eigenvalues for canonical form conditions
    IRI = ncon({gammaOut, lambdaOut, conj(gammaOut), conj(lambdaOut)},...
               {[1 3 2], [3 4], [1 5 2], [5 4]});
    norm2 = sqrt(IRI/D);
    gammaOut = gammaOut/norm2;
    
    mpsOut = {FORM, gammaOut, lambdaOut};
    normFactor = 1/(norm1*norm2);
end

function [mpsOut, AOut, normFactor] = AToCanonicalForm(mpsIn, TOL)
    % Transforms an mpsIn given in terms of A-tensors to canonical form 
    
    if nargin < 2
        TOL = 10;
    end
    
    % extract initial A tensor
    FORM = mpsIn{1};
    A    = mpsIn{2};
    
    if FORM ~= 1 && FORM ~= 2
        error('function called with wrong value of FORM = mpsIn{1}')
    end
    
    % Find bond dimension D
    D = size(A,1);

    % Construct transfer matrix E and reshape appropriately
    E = ncon({A,conj(A)},{[-1 -2 1],[-3 -4 1]});
    E = permute(E,[3,1,4,2]);
    E = reshape(E,[D^2,D^2]);

    % Finding the dominant eigenvectors of R and L, V_R and V_L respectively,
    % and the corresponding eigenvalues.
    [V_R,Eta_R]=eigs(E,1);
    [V_L,Eta_L]=eigs(E.',1);
    
    if abs(Eta_R - Eta_L) > 10^(-double(TOL))
        %E
        %['Eta_R - Eta_L = ',num2str(Eta_R - Eta_L)]
        error('Error: dominant left and right eigenvalues not equal');
    end

    % Divide by phase of first non-zero element to normalize
    idx_R = find(round(V_R, TOL),1);
    V_R = V_R/sign(V_R(idx_R)); % sign(z) = z./abs(z) for complex z
    
    idx_L = find(round(V_L, TOL),1);
    V_L = V_L/sign(V_L(idx_L));
    
    % Reshape V_R and V_L from D^2 element vectors to DxD matrices
    V_R=reshape(V_R,[D,D]).';
    V_L=reshape(V_L,[D,D]); % transposed in line with Orus2009 notation

    % raise error if either V_L or V_R aren't full rank
    VL_not_full_rank = rank(round(V_L, TOL)) < D;
    VR_not_full_rank = rank(round(V_R, TOL)) < D;
    
    if ( VL_not_full_rank || VR_not_full_rank )
        sprintf(' rank V_L = %d,',rank(round(V_L, TOL)));
        sprintf(' rank V_R = %d.',rank(round(V_R, TOL)));
        
        err.identifier = 'canonicalForm:VLVRnotFullRank';
        message1       = strcat('One or both of dominant right- and left-',...
                                ' eigenvectors of transfer matrix is',...
                                ' not full rank when reshaped to DxD matrix.');
        message2 = strcat(sprintf(' rank V_L = %d,',rank(round(V_L, TOL))), ...
                          sprintf(' rank V_R = %d.',rank(round(V_R, TOL))), ...
                          sprintf(' D = %d.', D), ...
                          sprintf(' V_R = '), mat2str(round(V_L,TOL)), ...
                          sprintf(' V_L = '), mat2str(round(V_R,TOL)));
        err.message = strcat(message1, message2);
        %error(err);
    end
    
    % Use squareDecomposition() to find X, Y
    X = squareDecomposition(V_R,TOL); % rightV = X*X'
    Y = squareDecomposition(V_L,TOL)';% leftV  = Y'*Y
    
    % Doing an SVD; Y*X = U*lambdaOut*V'
    [U,lambdaOut,V] = svd(Y*X);
    V=V'; % Redefine V, as svd() returns V and V' is needed
    
    % Perform final contractions dependent on whether right- or left-
    % canonical form is required
    switch mpsIn{1}
        case 1 % left canonical form
            AOut = ncon({lambdaOut,V,inv(X),A,inv(Y),U},...
                        {[-1 1],[1 2],[2 3],[3 4 -3],[4 5],[5 -2]});
            Lambda = ncon({lambdaOut, conj(lambdaOut)}, {[1 -1], [1 -2]});        
            
        case 2 % right canonical form
            AOut = ncon({V,inv(X),A,inv(Y),U,lambdaOut},...
                        {[-1 1],[1 2],[2 3 -3],[3 4],[4 5],[5,-2]});            
            Lambda = ncon({lambdaOut, conj(lambdaOut)}, {[-1 1], [-2 1]});
    end   
    
    % normalize so that tr(lambda.lambda')=1, and canonical form conditions
    % have unit eigenvalues
    % I. tr(lambda.lambda') = 1
    %Lambda = mpsOut{3};
    norm1 = trace(Lambda);
    Lambda = Lambda/norm1;
    
    % II. unit eigenvalues for canonical form conditions
    %IEI = ncon({AOut,conj(AOut)}, {[1 3 2], [1 3 2]});
    %norm2 = sqrt(IEI / D);
    norm2 = sqrt(Eta_R);
    AOut = AOut / norm2;
    
    mpsOut = {FORM, AOut, Lambda};
    normFactor = norm1*norm2;
end