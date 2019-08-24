function dA = inverseUsedTimestep( A, t, dt, H, KWARGS )
    %%INVERSEUSEDTIMESTEP Take either Euler, RK4, or symmetric timestep with
    %%inverse-using method
    
    f = @(a) findTimeDerivative(a, t, H, KWARGS);
    
    if KWARGS.RK4
        % 4th order Runge-Kutta
        k1 = f(A);
        k2 = f(A + 0.5*dt*k1);
        k3 = f(A + 0.5*dt*k2);
        k4 = f(A + dt*k3);
        
        dA = (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    elseif KWARGS.SYM_STEP
        dA = symmetricTimeStep(A, H, dt, size(A,1), KWARGS);
    else
        dA = dt*f(A); 
    end
end

%% Time-derivative and Symmetric step functions

function dAdt = findTimeDerivative(A, tiszero, H, KWARGS)
    [ ~, right, left, E ] = normalizeMPS(A);

    % Calculate V_L, the matrix of null vectors of lambda*A (L)
    if strcmp(KWARGS.Vmethod,'nullspace')
        sqrtleft = squareDecomposition(left, KWARGS.SQRT_TOL, 2);
        [~,D,d] = size(A);
        V_L = leftGaugeTangentBasis(sqrtleft, A, d, D, KWARGS.VL_TOL);            
    else
        if tiszero
            warning('Using unitary gauge, requires left canonical form');
            'This method is giving bad results!'
        end
        [A, right, U] = tangentBasis_unitaryGauge(A, right);
        U = permute(U, [1 2 4 3]); % U[i j rho sigma]
        V_L = reshape(U(:,:,2:end,:),[D D*(d-1) d]); % U[i (j,rho!=1), sigma]
    end


    % Compute C
    % NB for now considering one two-site operator h
    C = calculateC(A,H);

    % Calculate K
    K = calculateK(C,A,E,left,KWARGS.K_TOL);

    % Find inverses of l^(1/2), r^(1/2)
    lm12 = inv( squareDecomposition(left, KWARGS.SQRT_TOL, 2) );
    rm12 = inv( squareDecomposition(right,KWARGS.SQRT_TOL, 2) );

    % Calculate F
    F = calculateF(V_L,left,C,right,A,K,lm12,rm12,KWARGS.SQRT_TOL);

    % Calculate tangent space parameterization B(F)
    B_F = tangentVector(V_L,F,lm12,rm12,left,A,KWARGS.BF_TOL);

    dAdt = -1i*B_F; % -1i from Schrodinger equation
end


function dA = symmetricTimeStep(A, H, dt, D, KWARGS)
    % Take a symmetric timestep
    B     = @(a) findTimeDerivative( a, 0, H, KWARGS);
    cForm = @(a) canonicalForm({1,a,eye(D)},KWARGS.CFORM_TOL);
    nMax = 25;

    eta_threshold = 100*eps;
    eta = 10*eta_threshold;
    n = 0;

    % Iteratively look for midpoint
    A_midpoint_n = A;
    dAn = (dt/2) * B(A);

    while (eta > eta_threshold)
        % make nth guess
        A_midpoint_n = cForm(A_midpoint_n + dAn);  
        A_midpoint_n = A_midpoint_n{2};

        % calculate c, G, ln, rn
        [ c, G, l_n, r_n ] = midpointParams( A, A_midpoint_n, D);
        Ginv = inv(G);

        % Compute difference dAn between this and previous guess
        dAn = ncon({Ginv,A,G},{[-1 1],[1 2 -3],[2 -2]})/c - A_midpoint_n...
              + (dt/2) * B(A_midpoint_n);

        % Compute error measure
        eta = ncon({dAn, r_n, conj(dAn), l_n},...
                   {[1 2 3],[2 4],[5 4 3],[1 5]});
        % update n
        n = n + 1;

        if n==nMax
            error('unable to find midpoint');
        end            
    end 

    A_midpoint = cForm(A_midpoint_n + dAn);
    A_midpoint = A_midpoint{2};

    dA = A_midpoint + (dt/2) * B(A_midpoint) - A;
end


%% Subroutines for the above


function [ c, G, l_n, r_n ] = midpointParams( At, A_mp_n, D)
% Calculate various quantities needed for determining midpoint for timestep

    % find eigenstuff of E(A_midpoint_n) and E_{A_midpoint_n}^A
    E = ncon({A_mp_n, conj(A_mp_n)},...
                          {[-2 -4 1],[-1 -3 1]}); % E[i' i j' j]
    E = reshape(E,[D^2 D^2]);
    [l_n,~] = eigs(E.',1);
    [r_n,~] = eigs(E,1);

    EAA_midpoint_n = ncon({At, conj(A_mp_n)},...
                          {[-2 -4 1],[-1 -3 1]}); % E[i' i j' j]
    [l_nGinv, c] = eigs(reshape(EAA_midpoint_n,[D^2 D^2]).',1);
    
    % Divide by phase of first non-zero element to normalize
    idx_r = find(round(r_n, 8),1);
    r_n = r_n/sign(r_n(idx_r)); % sign(z) = z./abs(z) for complex z
    
    idx_l = find(round(l_n, 8),1);
    l_n = l_n/sign(l_n(idx_l));
    
    idx_l_nGinv = find(round(l_nGinv, 8),1);
    l_nGinv = l_nGinv/sign(l_nGinv(idx_l_nGinv));
    
    % reshape into matrices
    l_n = reshape(l_n, [D D]).'; 
    r_n = reshape(r_n, [D D]); 
    l_nGinv = reshape(l_nGinv, [D D]);
    
    % Work out G
    G = inv(inv(l_n).' * l_nGinv);
end


function C = calculateC(A,h)
    %% Calculate C, a [DxDxpxp] tensor
    % As shown in e.g. 'Summary' in Haegeman et al, PRB 88, 075133 (2013)
    %
    %% Inputs
    % A: [DxDxp] tensor
    %      MPS tensor
    % h: [pxpxpxp] tensor
    %      Hamiltonian tensor: h(i,j,k,l) = <i|<j|\hat{h}|k>|l>
    %% Diagram
    %          _____
    %  [-1] --|     |-- [-2]     [-1] -- A --- A -- [-2]
    %         |  C  |                    |_____|   
    %         |_____|         =         |       | 
    %          |   |                    |___H___|
    %        [-3] [-4]                   |     |   
    %%  
    
    % check dimensions
    p = size(A,3);
    if ~isequal(size(h), [p p p p])
        err.message = 'Hamiltonian tensor h has wrong dimension';
        err.identifier = 'calculateC:hWrongDimension';
        error(err);
    end
    
    % compute C
    C = ncon({A,A,h},{[-1 1 2], [1 -2 3], [-3 -4 2 3]});
    % (This entire fn could be replaced with one line, however having it as a 
    % separate function allows unit tests to be written for it)
end


function K=calculateK(C,A,E,left,K_TOL)
    %% Calculate K, a [DxD] matrix
    % As shown in e.g. 'Summary' in Haegeman et al, PRB 88, 075133 (2013)
    %
    %% Inputs
    % C: [DxDxpxp] tensor
    %      Two-site Hamiltonian tensor contracted with two A-tensors
    % A: [DxDxp] tensor
    %      MPS tensor normalized st dominant eigenvalue of transfer matrix = 1
    % E: [D^2xD^2] matrix
    %      Transfer matrix
    % left: [DxD] matrix
    %      Dominant left eigenvector of the transfer matrix E
    %
    %% Diagram:
    %
    %  /  |--      /  |--|     |--|     |--
    % |   |       |   |  |  C  |  |     |
    % | K |     = | l |  |_____|  | E_P |  
    % |   |       |   |   |   |   |     |
    %  \  |--      \  |-- A - A --|_____|--
    %% Another method (using iterative solver): - will be good for checking!
    % K = bicgstab(M,b) solves M*K=b for K
    % M = diag(ones(1,D)) - E + right*left
    % b = (ncon({C,A',A',left},{[1 -1 2 3],[4 5 2],[5 -2 3],[1 4]})) + ...
    % ((ncon({C,A',A',left,right},{[1 2 3 4],[5 6 3],[6 7 4],[1 5],[2 7]}))*left)
    %%
    
    if nargin==4
        K_TOL = 10;
    end
    
    % Find right and left eigenvectors and their associated eigenvalues
    D = sqrt(size(E,1));
    VR = zeros(D^2);
    VL = zeros(D^2);
    idx = 0; % degeneracies sometimes aren't normalized, so try up to 10 times
             % or until VR, VL are full-rank
    while ( rank(VR) < D^2  || rank(VL) < D^2 ) && idx < 10
        [VR,eigenvalues,VL] = eig(E); % cols of VL are e-vectors of E' (not E.')
        idx = idx + 1;
    end
    eigenvalues = diag(eigenvalues);
    maxEValue = max(eigenvalues);
    
    % If E has degeneracies eig() sometimes gives degenerate eigenvectors
    if ( rank(VR) < D^2 ) || ( rank(VL) < D^2 )
        err.identifier = 'TDVP:calculateK:degenerateEigenvectors';
        err.message =['Unable to find ',num2str(D^2),...
                      ' orthogonal eigenvectors after 10 iterations'];
        error(err)
    end
    
    % eig() gives us VL(:,i)'*VR(:,j!=i) = 0, however VL(:,i)'*VR(:,i) != 1,
    % so we need to normalize (but don't need to orthogonalize)
    
    for k=1:D^2
        norm = sqrt(VL(:,k)'*VR(:,k)); % norm is in general complex...
        VR(:,k) = VR(:,k)/norm;
        VL(:,k) = VL(:,k)/conj(norm);   % ...so divide by conjugate here
    end
    
    % check that sum_i |r_i)(l_i| = I
    overlap = @(i,j) VL(:,i)'*VR(:,j);
    iden = overlap(1:D^2, 1:D^2);
    if ~isequal(round(iden,K_TOL), eye(D^2))
        err = abs(max(max(iden - eye(D^2))));
        warning('sum_i |r_i)(l_i| = I + O(%.1d) != I',err);        
        if abs(err) > 1e-3
            [sprintf('iden = '),mat2str(iden),...
             sprintf('\nE = '),mat2str(E),...
             sprintf('\neigenvalues = '),mat2str(eigenvalues),...
             sprintf('\nVL = '),mat2str(VL),...
             sprintf('\nVR = '),mat2str(VR)]
        end
    end
                                             
    % E_P = sum_i [1/(1-lambda_i)]|ri)(li|
    E_P = zeros(D^2);
    for i=1:D^2
        if eigenvalues(i) ~= maxEValue
            coeff = 1 / (1 - eigenvalues(i));
            E_P=E_P + coeff * VR(:,i) * VL(:,i)'; 
            % NB we're using (v| = |v)' != |v).' here
        end
    end

    % reshape into DxDxDxD tensor
    E_P = permute(reshape(E_P,[D,D,D,D]), [2 4 1 3]);
               
    K = ncon({left, C, conj(A), conj(A), E_P},...
             {[1 2], [1 6 3 5], [2 4 3], [4 7 5], [6 -1 7 -2]});
end


function F = calculateF(V_L, left, C, right, A, K, lm12, rm12, TOL)
    %% Calculate F, a [(p-1)DxD] tensor
    % See e.g. 'Summary' in Haegeman et al, PRB 88, 075133 (2013)
    %
    %% Inputs
    % V_L:   [Dx(p-1)D] tensor
    %          columns of V_L form a basis over null space of A'*lambda tensor
    % left:  [DxD] matrix
    %          Dominant left eigenvector of the transfer matrix E
    % C:     [DxDxpxp] tensor
    %          Two-site Hamiltonian tensor contracted with two A-tensors
    % right: [DxD] matrix
    %          Dominant right eigenvector of the transfer matrix E
    % A:     [DxDxp] tensor
    %          MPS tensor
    % K:     [DxD] matrix
    %          tensor arising from sum of terms with B on RHS of Hamiltonian
    % lm12:  [DxD] matrix
    %          left^(-1/2) calculated previously
    % rm12:  [DxD] matrix
    %          right^(-1/2) calculated previously
    %%

    %% Compute x^(-1/2) and x^(1/2) for x = left, right
    lm12c = conj(lm12);
    rm12c = conj(rm12);

    l12 = squareDecomposition(left, TOL, 2);
    r12 = squareDecomposition(right, TOL, 2);
    
    %% Contract the three terms of F and sum
    
    % F1 = (l|H(AA,BA)|r)
    F1 = ncon({conj(V_L),l12,C,right,conj(A),rm12c},...
              {[2 -1 3], [1 2], [1 7 3 5], [7 6], [4 6 5], [-2 4]});
          
    % F2 = (l|H(AA,AB)|r)
    F2 = ncon({conj(V_L),lm12c.',conj(A),left,C,r12},...
              {[1 -1 2], [3 1], [5 3 4], [6 5], [6 7 4 2], [7 -2]});
 
    % F3 = sum_{k=0}^{inf} (l|H(AA,AA) E(A,A)^k E(A,B)|r)
    F3 = ncon({conj(V_L),lm12c.',K,A,r12},...
              {[1 -1 2], [3 1], [4 3], [4 5 2], [5 -2]});

    F = F1 + F2 + F3;
end


function BF = tangentVector(V_L,F,lm12,rm12,left,A,TOL)
    %% Calculate B_F, a [DxDxp] tensor
    % As shown in e.g. 'Summary' in Haegeman et al, PRB 88, 075133 (2013)
    %
    %%
    %    
    %    ________                           _______
    %  --| B(X) |--     --(left^(-1/2)).'-- | V_L | ~~ X -- (right^(-1/2)).' --
    %    |______|    =                      |_____|
    %       |                                  |
    %       |      
    %%
    
    BF = ncon({lm12.',V_L,F,rm12}, {[-1 1],[1 2 -3],[2 3],[3 -2]});
    
    % Adjust TOL by log10 of largest element of BF
    max_abs_el = abs(max(max(max(BF))));
    TOL = TOL - int64(ceil(log10(max_abs_el)));
    
    %% Check BF, V_L obey orthogonality as required
    D = size(left,1);
    should_be_zeros = ncon({conj(A),left,BF}, {[1 -1 2],[3 1],[3 -2 2]});
    if round(should_be_zeros, TOL) ~= zeros([D D])
        err.identifier = 'tangentVector:outputError';
        err.message = 'Output tensor is incorrect to within TOL';
        error(err)
    end

    should_be_zeros = ncon({conj(BF),left,A}, {[1 -1 2], [3 1], [3 -2 2]});
    if round(should_be_zeros, TOL) ~= zeros([D D])
        err.identifier = 'tangentVector:outputError';
        err.message = 'Output tensor is incorrect to within TOL';
        error(err)
    end    
end