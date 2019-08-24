function expectation = singleSiteExpectation( O, n, mps, Rn, Ln )
    %% Calculate expectation of operator O local to nth site 
    
    if nargin==3
        Rn_supplied = false;
        Ln_supplied = false;
    elseif nargin==4
        Rn_supplied = true;
        Ln_supplied = false;
    else
        Rn_supplied = true;
        Ln_supplied = true;
    end
    
    % 1. extract N, AN and check dims of O
    N = length(mps) - 1;
    
    An = mps{n+1}{1};
    if N==1 && size(An,1)==1, An = An.'; end
    
    dims = size(An);
    if n==1
        dn = dims(1);
    else
        dn = dims(end);
    end
    
    if ~isequal(size(O), [dn dn])
        errmsg = ['Operator on ',num2str(n),'th site has wrong dimensions:'];
        errmsg = [errmsg,newline,'size(O) = [',num2str(size(O)),'],',newline,...
                  'Local dimension = ',num2str(dn),'.'];
        error(errmsg)
    end
    
    % 2. calculate environments are required
    if ~Ln_supplied && ~Rn_supplied
        [Ln, Rn] = compute_env( mps, n, 'B' );
    elseif Ln_supplied && ~Rn_supplied
        [~, Rn] = compute_env( mps, n, 'R' );
    elseif ~Ln_supplied && Rn_supplied
        [Ln, ~] = compute_env( mps, n, 'L' );
    end
    
    % 3. Calculate <O>
    if n==1
        expectation = trace(An*Rn*An'*O);
    elseif n<N
        expectation = ncon({An,Ln,Rn,conj(An),O},{[1 3 2],[4 1],[3 6],[4 6 5],[5 2]});
    else
        expectation = trace(An*O.'*An'*Ln);
    end
    
    % 3. Raise warning if expectation has an imaginary part
    if abs(imag(expectation)) > 10*eps
        warning('Expecatation has an imaginary part')
    end
end

