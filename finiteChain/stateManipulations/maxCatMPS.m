function mpsOut = maxCatMPS( N, d )
    %%MAXCATMPS Outputs MPS of CAT state along cut with largest # Schmidt values
    %   This state has largest possible Schmidt values and may be useful
    %   when inverses are a problem in testing/debugging code.
    %
    %   mpsOut = maxCatMPS(N, d) 
    %       A max cat state over N sites each with local dimension d
    %   mpsOut = maxCatMPS(spinDimList) 
    %       A max cat state over a system with local dimensions given by the 
    %       array spinDimList = [d1,d2,...]
    %% Status: works for N,d input; doesn't for spinDimList input
    
    %% 1. Determine what kind of input we have and process
    if nargin == 2
        if ~isequal(size(N),[1 1]) || ~isequal(size(d),[1 1])
            error('Invalid input')
        end
        spinDimList = d*ones([1 N]);
    elseif nargin == 1
        spinDimList = N;
        if size(spinDimList,1)~=1
            error('Invalid input')
        end
        N = size(spinDimList,2);
    end        
    
    %% 2. If we have spinDimList input, find cut with max # Schmidt values
    if nargin == 1
        max_nSchmidtValues = 1;
        dimleft = 1;
        dimright = 1;
        for k = 1:N-1
            % loop over cuts
            dleft_k = prod(spinDimList(1:k));
            dright_k = prod(spinDimList(k+1:end));
            
            nSchmidtValues = min(dleft_k,dright_k);
            if nSchmidtValues > max_nSchmidtValues
                max_nSchmidtValues = nSchmidtValues;
                dimleft = dleft_k;
                dimright = dright_k;
            end
        end
    else
        n_middle = floor(N/2);
        dimleft = prod(spinDimList(1:n_middle));
        dimright = prod(spinDimList(n_middle+1:end));  
        max_nSchmidtValues = min(dimleft,dimright);
    end
    
    
    
    %% 3. Construct cat state over cut with max # Schmidt values
    
    max_entangled_1stcut = zeros([prod(spinDimList),1]); % state vector
    for k = 1:max_nSchmidtValues
        left_ket = zeros([1,dimleft]); left_ket(k) = 1;
        right_ket = zeros([1,dimright]); right_ket(k) = 1;
        max_entangled_1stcut = max_entangled_1stcut + ...
                               1/sqrt(max_nSchmidtValues)*kron(left_ket, right_ket).';
    end
    max_entangled_1stcut = reshape(max_entangled_1stcut,spinDimList);
    
    
    %% 4. Decompose into MPS   
    mpsOut = MPSdecomposition(max_entangled_1stcut);
    for k = 1:N-1 % print out Schmidt values along each cut
        sqrt(diag(mpsOut{k+1}{2})).';
    end
end

