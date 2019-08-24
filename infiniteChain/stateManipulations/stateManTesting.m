%% Main function to generate tests
% run by executing 'run(stateManTesting);' in the Command Window
function tests = stateManTesting
    tests = functiontests(localfunctions);
end

%% Test functions

function test_normalize_verifyOutput(testCase)
    % use AKLT mps (already in canonical form), multiply tensors by random
    % factors and attempt to renormalize
    TOL = 1e-10;
    
    dimRange = [2, 5];
    p = randi(dimRange);
    D = randi(dimRange);
    
    mps{2} = rand([D,D,p]) + rand([D,D,p])*1i;
    mps{3} = eye([D,D]);
        
    r1 = rand();  % random factors - could check we get these back
    r2 = rand();
    
    for FORM = 1:3
        mps{1} = FORM;
        mps = canonicalForm(mps); % this dep'nce isn't ideal, could use fixture
        
        notNormed = {FORM, mps{2}*r1, mps{3}*r2};
        
        normed = normalize(notNormed);
        
        % check inner product is 1
        expSolution = 1.;
        T = normed{2};
        L = normed{3};
        
        if FORM == 3
            L = L*L';
        end     
        
        % check sum of mod-square Schmidt values is one
        actSolution = trace(L);
        verifyEqual(testCase,actSolution,expSolution,'AbsTol',TOL);
        
        % NB could test more things eg. transfer matrix eigenvector, but I
        % think this is rigorous enough for now.
    end
end

function test_normalize_errMessage(testCase)
    % Test invalid FORM input throws the right error
    mpsIn = {6,[],[]};
    
    verifyError(testCase, @()normalize(mpsIn),...
            'normalize:InvalidFormInput');
end

function test_2siteDensityMatrix(testCase)
    % Test canonicalForm produces two-site density matrices that agree
    
    % Choose tolerances
    TOL = 1e-8;
    CFTOL = 8; % digits of tolerance for canonicalForm()
    
    % pick a p, D from a list of not-too-slow values
    dimRange = [2, 5];
    p = randi(dimRange);
    D = randi(dimRange);
    
    % generate random mps
    mpsTemp = mpsGenerator(p,D,3);
    [~,T,el] = mpsTemp{:}; % T, el are random and of the correct dimension
    
    % Find left-right eigenvectors of L, R transfer matrices
    [R, L] = constructELR({3,T,el});
    
    [leftV_R,eta_R] = eigs(R.',1);
    [rightV_L,eta_L] = eigs(L,1);
    
    leftV_R = reshape(leftV_R,[D,D]).';
    rightV_L = reshape(rightV_L,[D,D]).';
    
    if eta_R - eta_L > 1e-10
        error('got different right and left eigenvectors')
    end
    
    % normalize T so that rightV, leftV have eigenvalue 1 (and we can
    % therefore simplify the TN diagram)
    T = T/sqrt(eta_R);
    
    % find 1-site density matrix of this initial state
    numerator = ncon({T,el,T,rightV_L,leftV_R,conj(T),conj(el),conj(T)},...
                     {[1 2 -1], [2 3], [3 4 -2], [4 8],...
                      [1 5], [5 6 -3], [6 7], [7 8 -4]});
    denominator = ncon({T,rightV_L,leftV_R,conj(T)},...
                       {[1 2 5], [2 3], [1 4], [4 3 5]});
    
    dmInit = numerator/denominator;
    
    % find 2-site density matrix for each canonical form:
    % Vidal form
    mps = canonicalForm({3, T, el},CFTOL);
    [~,G,lv] = mps{:};
    dmVidal = ncon({lv,G,lv,G,lv,conj(lv),conj(G),conj(lv),conj(G),conj(lv)},...
                   {[1 2], [2 3 -1], [3 4], [4 5 -2], [5 6],...
                    [1 7], [7 8 -3], [8 9], [9 10 -4], [10 6]});
    
    % Left canonical form
    Lambda = ncon({el, conj(el)}, {[1 -1], [1 -2]});
    ALinit = ncon({el, T}, {[-1 1], [1 -2 -3]});
    mps = canonicalForm({1, ALinit, Lambda},CFTOL);
    [~,AL,LambdaL] = mps{:};
    dmLCanonical = ncon({AL,AL,LambdaL,conj(AL),conj(AL)},...
                        {[1 2 -1], [2 3 -2], [3 5], [1 4 -3], [4 5 -4]});
    
    % Right canonical form
    ARinit = ncon({T, el}, {[-1 2 -3], [2 -2]});
    mps = canonicalForm({2, ARinit, Lambda},CFTOL);
    [~,AR,LambdaR] = mps{:};
    dmRCanonical = ncon({AR,AR,conj(AR),conj(AR),LambdaR},...
                        {[1 2 -1], [2 3 -2], [4 5 -3], [5 3 -4], [1 4]});
    
    % Compare
    verifyEqual(testCase,dmInit,dmVidal,'AbsTol',TOL);
    verifyEqual(testCase,dmRCanonical,dmLCanonical,'AbsTol',TOL);
    verifyEqual(testCase,dmInit,dmRCanonical,'AbsTol',TOL);
end

function test_1siteDensityMatrix(testCase)
    % Test canonicalForm produces two-site density matrices that agree
    
    % Choose tolerance
    TOL = 1e-8;
    CFTOL = 8; % digits of tolerance in canonicalForm calls
    
    % pick a p, D from a list of not-too-slow values
    dimRange = [2, 5];
    p = randi(dimRange);
    D = randi(dimRange);
    
    % generate random mps
    mpsTemp = mpsGenerator(p,D,3);
    [~,T,el] = mpsTemp{:}; % T, el are random and of the correct dimension
    
    % Find left-right eigenvectors of L, R transfer matrices
    [R, L] = constructELR({3,T,el});
    
    [leftV_R,eta_R] = eigs(R.',1);
    [rightV_L,eta_L] = eigs(L,1);
    
    leftV_R = reshape(leftV_R,[D,D]).';
    rightV_L = reshape(rightV_L,[D,D]).';
    
    if eta_R - eta_L > 1e-10
        error('got different right and left eigenvectors')
    end
    
    % find 1-site density matrix of this initial state
    tensorList = {T,rightV_L,conj(T),leftV_R};
    numerator = ncon(tensorList,{[1 2 -1],[2 3],[4 3 -2], [1 4]});
    denominator = ncon(tensorList,{[1 2 5], [2 3], [4 3 5], [1 4]});
    
    dmInit = numerator/denominator;
    
    % find 1-site density matrix for each canonical form:
    % Vidal form
    mps = canonicalForm({3, T, el},CFTOL);
    [~,G,lambda] = mps{:};
    Lambdarght = lambda*lambda';
    Lambdaleft = lambda.'*conj(lambda);
    dmVidal = ncon({G,Lambdarght,conj(G),Lambdaleft},...
                           {[1 2 -1], [2 3], [4 3 -2], [1 4]});
    
    % Left canonical form
    Lambda = ncon({el, conj(el)}, {[1 -1], [1 -2]});
    ALinit = ncon({el, T}, {[-1 1], [1 -2 -3]});
    mps = canonicalForm({1, ALinit, Lambda},CFTOL);
    [~,AL,LambdaL] = mps{:};
    tensorList = {AL,LambdaL,conj(AL)};
    dmLCanonical = ncon(tensorList,{[1 2 -1], [2 3], [1 3 -2]});
    
    % Right canonical form
    ARinit = ncon({T, el}, {[-1 2 -3], [2 -2]});
    mps = canonicalForm({2, ARinit, Lambda},CFTOL);
    [~,AR,LambdaR] = mps{:};
    tensorList = {AR,conj(AR),LambdaR};
    dmRCanonical = ncon(tensorList,{[1 2 -1],[3 2 -2], [1 3]});
    
    % Compare
    verifyEqual(testCase,dmInit,dmVidal,'AbsTol',TOL);
    verifyEqual(testCase,dmRCanonical,dmLCanonical,'AbsTol',TOL);
    verifyEqual(testCase,dmInit,dmRCanonical,'AbsTol',TOL);
end

function test_SchmidtSpectrum(testCase)
    % Integration test to check that we get the same Schmidt spectrum of an
    % MPS whichever form we use to canonicalize
    
    % choose absolute tolerance
    TOL = 1e-10;
    CFTOL = 8; %digits of tolerance in canonicalForm calls
    
    % pick a p, D from a list of not-too-slow values
    dimRange = [2, 5];
    p = randi(dimRange);
    D = randi(dimRange);
    
    % generate random mps
    mpsTemp = mpsGenerator(p,D,3);
    [~,T,el] = mpsTemp{:}; % T, el are random and of the correct dimension
    
    % initialize empty cell-array-of-Schmidt-spectra
    schmidtSpectra = {0,0,0};
    
    for FORM = 1:3
        % put mps into form corresponding to the value of FORM
        if FORM == 3
            mps = {FORM, T, el};
        elseif FORM == 1
            Lambda = ncon({el, conj(el)}, {[1 -1], [1 -2]});
            ALinit = ncon({el, T}, {[-1 1], [1 -2 -3]});
            mps = {FORM, ALinit, Lambda};
        else
            Lambda = ncon({el, conj(el)}, {[1 -1], [1 -2]});
            ARinit = ncon({T, el}, {[-1 1 -3], [1 -2]});
            mps = {FORM, ARinit, Lambda};
        end
               
        % put into canonical form
        mps = canonicalForm(mps,CFTOL);
        
        % find Schmidt spectrum
        if FORM == 3
            schmidtSpectra{FORM} = diag(mps{3}); 
        else
            lambda = squareDecomposition(mps{3});
            schmidtSpectra{FORM} = diag(lambda); 
        end     
    end
    
    % check all three spectra are the same
    leftCanonical = schmidtSpectra{1};
    rightCanonical = schmidtSpectra{2};
    Vidal = schmidtSpectra{3};
    
    verifyEqual(testCase, Vidal, leftCanonical, 'AbsTol', TOL);
    verifyEqual(testCase, Vidal, rightCanonical, 'AbsTol', TOL);
    verifyEqual(testCase, leftCanonical, rightCanonical, 'AbsTol', TOL);
end

function test_Integrated1(testCase)
    % integration test to check that the functions
    %   mpsGenerator(),
    %   canonicalForm(), and
    %   verifyCanonicalForm(),
    % are working together properly.
    
    % index test cases by call-arrays [p,D,GLFORM,TOL1,TOL2]:
    %   TOL1 is used in GLToCanonicalForm and 
    %   TOL2 is used in verifyCanonicalForm
    check1 = {3, 2, 3, 10, 12}; % GLFORM always=3 for now since we have diff
    check2 = {2, 4, 3, 6, 9};   % functions for dealing with GL and A
    check3 = {5, 1, 3, 12, 12};   % (checked 1000 times in a for loop)
    
    checks = {check1, check2, check3};
    
    for k = 1:length(checks)
        [p,D,FORM,TOL1,TOL2] = checks{k}{:};
        state = mpsGenerator(p,D,FORM);
        cform = canonicalForm(state,TOL1);
        output = verifyCanonicalForm(cform,TOL2);
        verifyEqual(testCase,output,'TRUE');
    end
end

function test_canonicalForm_fixture_mpsAKLT(testCase)
    % test against 1d cluster state MPS
    clusterstate{2} = cat(3, [0 1;1 0], [0 -1i;1i 0], [1 0;0 -1]);
    clusterstate{3} = [1 0;0 1];
    
    TOL = 12;
    
    for k = 1:3 % check for each FORM
        clusterstate{1} = k;
        
        mps = canonicalForm(clusterstate);
        actSolution = round(mps{2},TOL);
        % get sign errors in the elements of G, which is ok, but need to be
        % corrected for
        actSolution(:,:,1) = actSolution(:,:,1)/sign(real(actSolution(1,2,1)));
        actSolution(:,:,2) = actSolution(:,:,2)/sign(imag(actSolution(2,1,2)));
        actSolution(:,:,3) = actSolution(:,:,3)/sign(real(actSolution(1,1,3)));

        if k ~= 0
            % the extra tensor lambda in A compared to Gammacontributes a 
            % factor 1/sqrt(2) we need to correct for
            actSolution = sqrt(2) * actSolution;
        end
        
        expSolution = sqrt(2)*clusterstate{2};
        expSolution = round(expSolution,TOL); 

        verifyEqual(testCase,actSolution,expSolution,'RelTol',TOL);
    end
end
    
function test_mpsGenerator_AFORM(testCase)
    % test that a cell-array of the correct dimensions, containing arrays
    % of the correct dimensions, is produced when calling with FORM=1,2

    for FORM = 1:2
        state = mpsGenerator(6,5,FORM);
        A = state{2};

        actstateSize = size(state);
        expstateSize = [1 3];
        verifyEqual(testCase,actstateSize,expstateSize);

        actTensorDimns = size(A);
        expTensorDimns = [5 5 6];
        verifyEqual(testCase,actTensorDimns,expTensorDimns);
    end
end

function test_mpsAKLT(testCase)
    % test we're getting the right output
    state = mpsAKLT(0);
    actSolution = state{2};
    expSolution = cat(3, [0 1;1 0], [0 -1i;1i 0], [1 0;0 -1]);
    
    verifyEqual(testCase,actSolution,expSolution)
end

function test_verifyCanonicalForm(testCase)
    % test that two example cases are correctly verified
    G = cat(3,[0 0;0 1.92066023856587],[1.79551187099548 0;0 0],...
              [1.79551187099548 0;0 1.92066023856587]);
    l = [0.393819051051435 0;0 0.368158181748488];
    mps = {3,G,l};
    
    verifyEqual(testCase,verifyCanonicalForm(mps),'TRUE');
    
    G = cat(3, [0 1;1 0], [0 -1i;1i 0], [0 0;0 -1]);
    l = [1 0;0 1];
    mps2 = {3,G,l};
    
    verifyEqual(testCase,verifyCanonicalForm(mps2),'FALSE');
end

function test_squareDecomposition_type1(testCase)
    % test output works for an example
    M = [4 0 3;0 4 1;3 1 3]; % a Hermitian, posdef matrix
    X = squareDecomposition(M);
    
    expOutput = round(M,12);
    actOutput = round(X*X',12);
    
    verifyEqual(testCase,actOutput,expOutput);
end

function test_squareDecomposition_type2(testCase)
    % test output works for an example
    M = [4 0 3;0 4 1;3 1 3]; % a Hermitian, posdef matrix
    X = squareDecomposition(M,10,2);
    
    expOutput = round(M,12);
    actOutput = round(X*X,12);
    
    verifyEqual(testCase,actOutput,expOutput);
end
    
function [R, L] = constructELR(mpsIn)
    %  constructs from MPS the appropriate of E, [R,L] with manual
    %  reshaping
    FORM = mpsIn{1};
    
    if FORM == 2
        FORM = 1;
    end
    
    [D,~,~] = size(mpsIn{2});
    
    switch FORM
        case 3
            [~,T,el] = mpsIn{:};
            % Find R and L matrices as rank-4 tensors
            tensors = {T, el, conj(T), conj(el)};
            legLinksRight = {[-1 2 1], [2 -2], [-3 3 1], [3 -4]};
            legLinksLeft  = {[2 -2 1], [-1 2], [3 -4 1], [-3 3]};
            R = ncon(tensors, legLinksRight); % order of indices in R, L is
            L = ncon(tensors, legLinksLeft);  % alpha, beta, alpha', beta'

            % Reshape R and L vectors explicitly by looping over each new index 
            newR = zeros([D^2,D^2]);
            newL = zeros([D^2,D^2]);

            % Loop over the new indices
            for aa=1:D^2       % aa ie. (alpha alpha')
                for bb=1:D^2   % bb ie. (beta beta')
                    a  = floor(aa./(D + 1e-5)) + 1;
                    ap = aa - (a - 1)*D;
                    b  = floor(bb./(D + 1e-5)) + 1;
                    bp = bb - (b - 1)*D;

                    % assign elements of newR, newL correspondingly
                    newR(aa, bb) = R(a, b, ap, bp);
                    newL(aa, bb) = L(a, b, ap, bp);
                end
            end

            % Set R and L equal to their reshaped counterparts
            R = newR;
            L = newL;
        case 1
            [~,A,~] = mpsIn{:};
            newE = zeros([D^2, D^2]);
            E = ncon({A,conj(A)},{[-1 -2 1],[-3 -4 1]});
            for aa=1:D^2       % aa ie. (alpha alpha')
                for bb=1:D^2   % bb ie. (beta beta')
                    a  = floor(aa./(D + 1e-5)) + 1;
                    ap = aa - (a - 1)*D;
                    b  = floor(bb./(D + 1e-5)) + 1;
                    bp = bb - (b - 1)*D;

                    % assign elements of newE correspondingly
                    newE(aa, bb) = E(a, b, ap, bp);
                end
            end

            E = newE;
            R = E;
            L = E;
    end
end