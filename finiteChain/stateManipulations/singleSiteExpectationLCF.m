function expectation = singleSiteExpectationLCF( O, n, mps_LCF )
    %% Calculate expectation of operator O local to nth site with mps_LCF
    
    % 1. extract variables and check dims of O
    N = length(mps_LCF) - 1;
    
    An = mps_LCF{n+1}{1};
    if N==1 && size(An,1)==1, An = An.'; end
    Rn = mps_LCF{n+1}{2};
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
    
    % 2. Calculate <O>
    if n==1
        expectation = trace(An*Rn*An'*O);
    elseif n<N
        expectation = ncon({An,Rn,conj(An),O},{[1 2 5],[2 3],[1 3 4],[4 5]});
    else
        expectation = trace(An*O.'*An');
    end
    
    % 3. Raise warning if expectation has an imaginary part
    if abs(imag(expectation)) > 10*eps
        warning('Expecatation has an imaginary part')
    end
end

