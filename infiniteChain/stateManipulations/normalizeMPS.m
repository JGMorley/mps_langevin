function [ A, r, l, E ] = normalizeMPS( mpsIn )
    %% Normalizes A and returns this with normalized dominant eigenvectors r, l
    %
    %% Inputs
    % mpsIn: cell-array or tensor A
    %%
    
    % Find dominant eigenvalue of transfer matrix, eta
    if iscell(mpsIn)
        A = mpsIn{2};
    else
        A = mpsIn;
    end
        
    if size(size(A),2) ~= 3
        err.identifier = 'normalizeMPS:incorrectRankA';
        err.message = 'Input MPS tensor A is not of rank 3';
        error(err)
    end
    
    D = size(A,1);
    E = ncon({A,conj(A)},{[-2 -4 1],[-1 -3 1]}); % E[i' i j' j]
    E = reshape(E,[D^2,D^2]); % E[(i i') (j j')]
    
    % Check we have a unique dominant eigenvalue if D>1
    if D > 1
        [~, eta12] = eigs(E,2); % find first two eigenvalues/vectors
        diffAbsValues = abs(abs(eta12(1,1)) - abs(eta12(2,2)));
        if diffAbsValues < eps
            warning('More than one dominant eigenvalue');
        end
    end
    % Find dominant left, right eigenvectors
    [r,eta] = eigs(E,1);
    [l,~] = eigs(E.',1); % Tranpose so we can contract lE = l
    
    % Reshape and normalize
    A = A / sqrt(eta);
    E = E / eta;
    
    overlap = l.'*r;
    r = reshape(r,[D,D]).' / sqrt(overlap); 
    l = reshape(l,[D,D]).' / sqrt(overlap); % so lE = l and tr(l.'*r) = 1
    %l = reshape(l,[D,D]).' / sqrt(overlap); % so l.'E = l.' and tr(l*r) = 1
    
    % Correct for phase - l, r are posdef up to an overall phase
    carr = {r, l};
    for k = 1:2
        x = carr{k};   
        maxAbs = max(max(abs(x))); % max() works column-wise
        idx = find(abs(x)==maxAbs(1));
        phase = sign(x(idx(1))); % phase of largest magnitude element     
        carr{k} = x / phase;            
    end
    [r,l] = carr{:};
    
    % If we're in canonical form, renormalize r and l so that the
    % appropriate environment is the identity, not just up to normalization
    if strcmp(verifyCanonicalForm({1,A,r}),'TRUE') 
        % LCF
        norm_factor = l(1);
        l = l/norm_factor;
        r = r*norm_factor;
    elseif strcmp(verifyCanonicalForm({2,A,l}),'TRUE')
        % RCF
        norm_factor = r(1);
        l = l*norm_factor;
        r = r/norm_factor;
    end
end

%% check with eg.
% A = rand([3,3,2]) + rand([3,3,2])*1i;
% [A,R,L] = normalizeMPS(A);
% trace(L.'*R); % should equal 1
% E = ncon({A,conj(A)},{[-2 -4 1],[-1 -3 1]});
% E = reshape(E,[D^2,D^2]);
% eigs(E,1); % should equal 1
% E = permute(reshape(E, [D D D D]), [2 1 4 3]);
% ncon({L, E}, {[1 2], [1 2 -1 -2]}) - L % should get all zeros to eps
% ncon({E, R}, {[-1 -2 1 2], [1 2]}) - R % should get all zeros to eps
